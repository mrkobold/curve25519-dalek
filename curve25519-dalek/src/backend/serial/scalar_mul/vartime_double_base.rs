// -*- mode: rust; -*-
//
// This file is part of curve25519-dalek.
// Copyright (c) 2016-2021 isis lovecruft
// Copyright (c) 2016-2019 Henry de Valence
// See LICENSE for licensing information.
//
// Authors:
// - isis agora lovecruft <isis@patternsinthevoid.net>
// - Henry de Valence <hdevalence@hdevalence.ca>
#![allow(non_snake_case)]

use core::cmp::Ordering;
extern crate std;
use std::string::ToString;

use crate::backend::serial::curve_models::{ProjectiveNielsPoint, ProjectivePoint};
use crate::backend::serial::u64::field::FieldElement51;
use crate::constants;
use crate::edwards::EdwardsPoint;
use crate::scalar::Scalar;
use crate::traits::Identity;
use crate::window::NafLookupTable5;

/// Compute \\(aA + bB\\) in variable time, where \\(B\\) is the Ed25519 basepoint.
pub fn mul(a: &Scalar, A: &EdwardsPoint, b: &Scalar) -> EdwardsPoint {
    let a_naf = a.non_adjacent_form(5);

    #[cfg(feature = "precomputed-tables")]
    let b_naf = b.non_adjacent_form(8);
    #[cfg(not(feature = "precomputed-tables"))]
    let b_naf = b.non_adjacent_form(5);

    // Find starting index
    let mut i: usize = 255;
    for j in (0..256).rev() {
        i = j;
        if a_naf[i] != 0 || b_naf[i] != 0 {
            break;
        }
    }

    let table_A = NafLookupTable5::<ProjectiveNielsPoint>::from(A);
    #[cfg(feature = "precomputed-tables")]
    let table_B = &constants::AFFINE_ODD_MULTIPLES_OF_BASEPOINT;
    #[cfg(not(feature = "precomputed-tables"))]
    let table_B =
        &NafLookupTable5::<ProjectiveNielsPoint>::from(&constants::ED25519_BASEPOINT_POINT);

    let mut r = ProjectivePoint::identity();
    loop {
        let mut t = r.double();

        match a_naf[i].cmp(&0) {
            Ordering::Greater => t = &t.as_extended() + &table_A.select(a_naf[i] as usize),
            Ordering::Less => t = &t.as_extended() - &table_A.select(-a_naf[i] as usize),
            Ordering::Equal => {}
        }

        match b_naf[i].cmp(&0) {
            Ordering::Greater => t = &t.as_extended() + &table_B.select(b_naf[i] as usize),
            Ordering::Less => t = &t.as_extended() - &table_B.select(-b_naf[i] as usize),
            Ordering::Equal => {}
        }

        r = t.as_projective();

        if i == 0 {
            break;
        }
        i -= 1;
    }

    r.as_extended()
}

pub fn deserialize_r_from_backup(projective_point_bu: [u64; 15]) -> ProjectivePoint {
    let mut array: [u64; 5] = [0u64; 5];
    array.clone_from_slice(&projective_point_bu[0..5]);
    let x_field_element = FieldElement51(array);

    let mut array: [u64; 5] = [0u64; 5];
    array.clone_from_slice(&projective_point_bu[5..10]);
    let y_field_element = FieldElement51(array);

    let mut array: [u64; 5] = [0u64; 5];
    array.clone_from_slice(&projective_point_bu[10..15]);
    let z_field_element = FieldElement51(array);

    ProjectivePoint {
        X: x_field_element,
        Y: y_field_element,
        Z: z_field_element,
    }
}

/// Compute \\(aA + bB\\) in variable time, where \\(B\\) is the Ed25519 basepoint.
pub fn kobold_mul<F: Fn(usize, [u64; 15]) -> ()>(
    a: &Scalar,
    A: &EdwardsPoint,
    b: &Scalar,
    update_kobold_account_handle: F,
    i_bu: usize,
    projective_point_bu: [u64; 15],
) -> (EdwardsPoint, u8) {
    let a_naf = a.non_adjacent_form(5);

    #[cfg(feature = "precomputed-tables")]
    let b_naf = b.non_adjacent_form(8);
    #[cfg(not(feature = "precomputed-tables"))]
    let b_naf = b.non_adjacent_form(5);

    let mut i: usize = 255;
    let mut r;
    if i_bu == 300 {
        // Find starting index
        for j in (0..256).rev() {
            i = j;
            if a_naf[i] != 0 || b_naf[i] != 0 {
                break;
            }
        }
        r = ProjectivePoint::identity();
    } else {
        i = i_bu;
        r = deserialize_r_from_backup(projective_point_bu);
    }

    let table_A = NafLookupTable5::<ProjectiveNielsPoint>::from(A);
    #[cfg(feature = "precomputed-tables")]
    let table_B = &constants::AFFINE_ODD_MULTIPLES_OF_BASEPOINT;
    #[cfg(not(feature = "precomputed-tables"))]
    let table_B =
        &NafLookupTable5::<ProjectiveNielsPoint>::from(&constants::ED25519_BASEPOINT_POINT);

    let mut j = 0;
    loop {
        let mut t = r.double();

        match a_naf[i].cmp(&0) {
            Ordering::Greater => t = &t.as_extended() + &table_A.select(a_naf[i] as usize),
            Ordering::Less => t = &t.as_extended() - &table_A.select(-a_naf[i] as usize),
            Ordering::Equal => {}
        }

        match b_naf[i].cmp(&0) {
            Ordering::Greater => t = &t.as_extended() + &table_B.select(b_naf[i] as usize),
            Ordering::Less => t = &t.as_extended() - &table_B.select(-b_naf[i] as usize),
            Ordering::Equal => {}
        }

        r = t.as_projective();

        if i == 0 {
            break;
        }
        if j == 30 {
            i -= 1;
            break;
        }
        i -= 1;
        j += 1;
    }

    // i == 0 and the algo already performed the repetition for i = 0
    if i == 0 && j < 30 {
        return (r.as_extended(), 2);
    } else {
        // save progress
        let x_arr = r.X.0;
        let y_arr = r.Y.0;
        let z_arr = r.Z.0;
        let concatenated = [x_arr, y_arr, z_arr].concat();
        let projective_point_progress: [u64; 15] = concatenated.try_into().unwrap(); // TODO: verify correctness of this
        update_kobold_account_handle(i, projective_point_progress);
        return (EdwardsPoint::default(), 1);
    }
}
