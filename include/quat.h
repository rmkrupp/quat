/* File: include/quat.h
 * Part of quat <github.com/rmkrupp/quat>
 *
 * Copyright (C) 2024 Noah Santer <n.ed.santer@gmail.com>
 * Copyright (C) 2024 Rebecca Krupp <beka.krupp@gmail.com>
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
#ifndef QUAT_H
#define QUAT_H

/*
 * A SIMPLE QUATERNION & MATRIX LIBRARY
 */

/* Notes
 *
 * for all these functions, it is perfectly fine for the output pointer to
 * be the same as any of the input pointers of the same type
 *
 * all matrix arguments are 4x4 (16-element) matrices
 */

typedef double quaternion_element;
typedef float matrix_element;

struct matrix {
    matrix_element matrix[16];
};

struct quaternion {
    quaternion_element x,
                       y,
                       z,
                       w;
};

void quaterion_identity(
        struct quaternion * out
    ) [[gnu::nonnull(1)]];

void quaternion_mul(
        struct quaternion * out,
        const struct quaternion * a,
        const struct quaternion * b
    ) [[gnu::nonnull(1, 2, 3)]];

void quaternion_normalize(
        struct quaternion * out,
        const struct quaternion * a
    ) [[gnu::nonnull(1, 2)]];

void quaternion_conjugate(
        struct quaternion * out,
        const struct quaternion * a
    ) [[gnu::nonnull(1, 2)]];

void quaternion_to_matrix(
        struct matrix * matrix,
        const struct quaternion * a
    ) [[gnu::nonnull(1, 2)]];

void quaternion_from_axis_angle(
        struct quaternion * out,
        double ax,
        double ay,
        double az,
        double theta
    ) [[gnu::nonnull(1)]];

void matrix_identity(
        struct matrix * out
    ) [[gnu::nonnull(1)]];

void matrix_translation(
        struct matrix * out,
        matrix_element x,
        matrix_element y,
        matrix_element z
    ) [[gnu::nonnull(1)]];

void matrix_perspective(
        struct matrix * out,
        matrix_element z_far,
        matrix_element z_near,
        matrix_element f_scale,
        matrix_element aspect
    ) [[gnu::nonnull(1)]];

void matrix_multiply(
        struct matrix * out,
        const struct matrix * a,
        const struct matrix * b
    ) [[gnu::nonnull(1, 2, 3)]];

#endif /* HASH_H */
