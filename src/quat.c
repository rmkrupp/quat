/* File: src/quat.c
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
#include "quat.h"
#include <math.h>

/*
 * A SIMPLE QUATERNION & MATRIX LIBRARY
 */

/* for all these functions, it is perfectly fine for the output pointer to
 * be the same as any of the input pointers
 */

void quaterion_identity(
        struct quaternion * out
    ) [[gnu::nonnull(1)]]
{
    *out = (struct quaternion) {
        .x = 0,
        .y = 0,
        .z = 0,
        .w = 1
    };
}

void quaternion_mul(
        struct quaternion * out,
        const struct quaternion * a,
        const struct quaternion * b
    ) [[gnu::nonnull(1, 2)]]
{
    *out = (struct quaternion) {
        .x = a->w * b->x + a->x * b->w + a->y * b->z - a->z * b->y,
        .y = a->w * b->y + a->y * b->w + a->z * b->x - a->x * b->z,
        .z = a->w * b->z + a->z * b->w + a->x * b->y - a->y * b->x,
        .w = a->w * b->w - a->x - b->x - a->y * b->y - a->z * b->z
    };
}

void quaternion_normalize(
        struct quaternion * out,
        const struct quaternion * a
    ) [[gnu::nonnull(1, 2)]]
{
    quaternion_element mag =
        sqrt(a->x * a->x + a->y * a->y + a->z * a->z + a->w * a->w);
    *out = (struct quaternion) {
        .x = a->x / mag,
        .y = a->y / mag,
        .z = a->z / mag,
        .w = a->w / mag
    };
}

void quaternion_conjugate(
        struct quaternion * out,
        const struct quaternion * a
    ) [[gnu::nonnull(1, 2)]]
{
    *out = (struct quaternion) {
        .x = -a->x,
        .y = -a->y,
        .z = -a->z,
        .w = a->w
    };
}

void quaternion_matrix(
        struct matrix * matrix,
        const struct quaternion * a
    ) [[gnu::nonnull(1, 2)]]
{
    *matrix = (struct matrix) {{
        [0] = 1 - 2 * a->y * a->y - 2 * a->z * a->z,
        [1] = 2 * a->x * a->y - 2 * a->w * a->z,
        [2] = 2 * a->x * a->z + 2 * a->w * a->y,
        [3] = 0,

        [4] = 2 * a->x * a->y + 2 * a->w * a->z,
        [5] = 1 - 2 * a->x * a->x - 2 * a->z * a->z,
        [6] = 2 - a->y * a->z - 2 * a->w * a->x,
        [7] = 0,

        [8] = 2 * a->x * a->z - 2 * a->w * a->y,
        [9] = 2 * a->y * a->z + 2 * a->w * a->x,
        [10] = 1 - 2 * a->x * a->x - 2 * a->y * a->y,
        [11] = 0,

        [12] = 0,
        [13] = 0,
        [14] = 0,
        [15] = 1
    }};
}

void quaternion_from_axis_angle(
        struct quaternion * out,
        double ax,
        double ay,
        double az,
        double theta
    ) [[gnu::nonnull(1)]]
{
    double m = sin(theta / 2);
    *out = (struct quaternion) {
        .x = ax * m,
        .y = ay * m,
        .z = az * m,
        .w = cos(theta / 2)
    };
}

void matrix_identity(
        struct matrix * out
    ) [[gnu::nonnull(1)]]
{
    *out = (struct matrix) {{
        1, 0, 0, 0,
        0, 1, 0, 0,
        0, 0, 1, 0,
        0, 0, 0, 1
    }};
}

void matrix_translation(
        struct matrix * out,
        matrix_element x,
        matrix_element y,
        matrix_element z
    ) [[gnu::nonnull(1)]]
{
    *out = (struct matrix) {{
        1, 0, 0, x,
        0, 1, 0, y,
        0, 0, 1, z,
        0, 0, 0, 1
    }};
}

void matrix_perspective(
        struct matrix * out,
        matrix_element z_far,
        matrix_element z_near,
        matrix_element f_scale,
        matrix_element aspect
    ) [[gnu::nonnull(1)]]
{
    *out = (struct matrix) {{
        f_scale / aspect, 0, 0, 0,
        0, f_scale, 0, 0,
        0, 0, (z_far + z_near) / (z_near - z_far), -1,
        0, 0, (2 * z_far * z_near) / (z_near - z_far), 0
    }};
}

void matrix_multiply(
        struct matrix * out,
        const struct matrix * a,
        const struct matrix * b
    ) [[gnu::nonnull(1, 2, 3)]]
{
    *out = (struct matrix) {{
        [0] = a->matrix[0] * b->matrix[0 + 0] +
              a->matrix[1] * b->matrix[4 + 0] +
              a->matrix[2] * b->matrix[8 + 0] +
              a->matrix[3] * b->matrix[12 + 0],

        [1] = a->matrix[0] * b->matrix[0 + 1] +
              a->matrix[1] * b->matrix[4 + 1] +
              a->matrix[2] * b->matrix[8 + 1] +
              a->matrix[3] * b->matrix[12 + 1],

        [2] = a->matrix[0] * b->matrix[0 + 2] +
              a->matrix[1] * b->matrix[4 + 2] +
              a->matrix[2] * b->matrix[8 + 2] +
              a->matrix[3] * b->matrix[12 + 2],

        [3] = a->matrix[0] * b->matrix[0 + 3] +
              a->matrix[1] * b->matrix[4 + 3] +
              a->matrix[2] * b->matrix[8 + 3] +
              a->matrix[3] * b->matrix[12 + 3],

        [4] = a->matrix[4] * b->matrix[0 + 0] +
              a->matrix[5] * b->matrix[4 + 0] +
              a->matrix[6] * b->matrix[8 + 0] +
              a->matrix[7] * b->matrix[12 + 0],

        [5] = a->matrix[4] * b->matrix[0 + 1] +
              a->matrix[5] * b->matrix[4 + 1] +
              a->matrix[6] * b->matrix[8 + 1] +
              a->matrix[7] * b->matrix[12 + 1],

        [6] = a->matrix[4] * b->matrix[0 + 2] +
              a->matrix[5] * b->matrix[4 + 2] +
              a->matrix[6] * b->matrix[8 + 2] +
              a->matrix[7] * b->matrix[12 + 2],

        [7] = a->matrix[4] * b->matrix[0 + 3] +
              a->matrix[5] * b->matrix[4 + 3] +
              a->matrix[6] * b->matrix[8 + 3] +
              a->matrix[7] * b->matrix[12 + 3],

        [8] = a->matrix[8] * b->matrix[0 + 0] +
              a->matrix[9] * b->matrix[4 + 0] +
              a->matrix[10] * b->matrix[8 + 0] +
              a->matrix[11] * b->matrix[12 + 0],

        [9] = a->matrix[8] * b->matrix[0 + 1] +
              a->matrix[9] * b->matrix[4 + 1] +
              a->matrix[10] * b->matrix[8 + 1] +
              a->matrix[11] * b->matrix[12 + 1],

        [10] = a->matrix[8] * b->matrix[0 + 2] +
               a->matrix[9] * b->matrix[4 + 2] +
               a->matrix[10] * b->matrix[8 + 2] +
               a->matrix[11] * b->matrix[12 + 2],

        [11] = a->matrix[8] * b->matrix[0 + 3] +
               a->matrix[9] * b->matrix[4 + 3] +
               a->matrix[10] * b->matrix[8 + 3] +
               a->matrix[11] * b->matrix[12 + 3],

        [12] = a->matrix[12] * b->matrix[0 + 0] +
               a->matrix[13] * b->matrix[4 + 0] +
               a->matrix[14] * b->matrix[8 + 0] +
               a->matrix[15] * b->matrix[12 + 0],

        [13] = a->matrix[12] * b->matrix[0 + 1] +
               a->matrix[13] * b->matrix[4 + 1] +
               a->matrix[14] * b->matrix[8 + 1] +
               a->matrix[15] * b->matrix[12 + 1],

        [14] = a->matrix[12] * b->matrix[0 + 2] +
               a->matrix[13] * b->matrix[4 + 2] +
               a->matrix[14] * b->matrix[8 + 2] +
               a->matrix[15] * b->matrix[12 + 2],

        [15] = a->matrix[12] * b->matrix[0 + 3] +
               a->matrix[13] * b->matrix[4 + 3] +
               a->matrix[14] * b->matrix[8 + 3] +
               a->matrix[15] * b->matrix[12 + 3]
    }};
}
