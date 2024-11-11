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
 */

/* quaternions are four-vectors made of this element type */
typedef double quaternion_element;

/* matrices are sixteen-vectors (4x4) made of this element type */
typedef float matrix_element;

/* a 4x4 matrix */
struct matrix {
    matrix_element matrix[16];
};

/* a quaternion */
struct quaternion {
    quaternion_element x,
                       y,
                       z,
                       w;
};

/* the magnitude (length) of the quaternion; if 1, it is a unit quaternion */
quaternion_element quaternion_magnitude(
        const struct quaternion * a
    ) [[gnu::nonnull(1)]];

/* set out to the identity quaternion <0 0 0 1> */
void quaterion_identity(
        struct quaternion * out
    ) [[gnu::nonnull(1)]];

/* set out to the product of a and b */
void quaternion_mul(
        struct quaternion * out,
        const struct quaternion * a,
        const struct quaternion * b
    ) [[gnu::nonnull(1, 2, 3)]];

/* set out to a, normalized: out is a with all components divided by m such
 * that the sum of the squares of components of out is 1
 */
void quaternion_normalize(
        struct quaternion * out,
        const struct quaternion * a
    ) [[gnu::nonnull(1, 2)]];

/* set out to the conjugate of a. if a is a unit quaternion, the conjugate
 * is its inverse
 */
void quaternion_conjugate(
        struct quaternion * out,
        const struct quaternion * a
    ) [[gnu::nonnull(1, 2)]];

/* set matrix to the matrix form of a */
void quaternion_to_matrix(
        struct matrix * matrix,
        const struct quaternion * a
    ) [[gnu::nonnull(1, 2)]];

/* set out to the quaternion calculated from rotation theta around axis
 * (ax, ay, az)
 */
void quaternion_from_axis_angle(
        struct quaternion * out,
        double ax,
        double ay,
        double az,
        double theta
    ) [[gnu::nonnull(1)]];

/* set out to the identiy matrix */
void matrix_identity(
        struct matrix * out
    ) [[gnu::nonnull(1)]];

/* set out to a translation matrix for vector <x y z> */
void matrix_translation(
        struct matrix * out,
        matrix_element x,
        matrix_element y,
        matrix_element z
    ) [[gnu::nonnull(1)]];

/* set out to a perspective matrix with these properties */
void matrix_perspective(
        struct matrix * out,
        matrix_element z_far,
        matrix_element z_near,
        matrix_element f_scale,
        matrix_element aspect
    ) [[gnu::nonnull(1)]];

/* set out to the product of a and b */
void matrix_multiply(
        struct matrix * out,
        const struct matrix * a,
        const struct matrix * b
    ) [[gnu::nonnull(1, 2, 3)]];

#endif /* HASH_H */
