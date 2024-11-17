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

/* the magnitude (length) of the quaternion; if 1, it is a unit quaternion */
quaternion_element quaternion_magnitude(
        const struct quaternion * a
    ) [[gnu::nonnull(1)]]
{
    return sqrt(a->x * a->x + a->y * a->y + a->z * a->z + a->w * a->w);
}

/* set out to the identity quaternion <0 0 0 1> */
void quaternion_identity(
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

/* set out to the product of a and b
 *
 * this is equivalent to:
 *
 * vec3 s1v2 = a.w * <b.x b.y b.z>
 * vec3 s2v1 = b.w * <a.x a.y a.z>
 * vec3 cross = <a.x a.y a.z> x <b.x b.y b.z>
 * out = {
 *     .x = cross.x + a.w * b.x + b.w * a.x,
 *     .y = cross.y + a.w * b.y + b.w * a.y,
 *     .z = cross.z + a.w * b.z + b.w * a.z,
 *     .w = a.w * b.w - a.x * b.x - a.y * b.y - a.z * b.z
 * }
 *
 */
void quaternion_multiply(
        struct quaternion * out,
        const struct quaternion * a,
        const struct quaternion * b
    ) [[gnu::nonnull(1, 2)]]
{
    *out = (struct quaternion) {
        .x = a->w * b->x + b->w * a->x + a->y * b->z - a->z * b->y,
        .y = a->w * b->y + b->w * a->y - a->x * b->z + a->z * b->x,
        .z = a->w * b->z + b->w * a->z + a->x * b->y - a->y * b->x,
        .w = a->w * b->w - a->x * b->x - a->y * b->y - a->z * b->z
    };
}


/* set out to a, normalized: out is a with all components divided by m such
 * that the sum of the squares of components of out is 1
 */
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

/* set out to the conjugate of a. if a is a unit quaternion, the conjugate
 * is its inverse
 */
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

/* set matrix to the matrix form of a */
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
        [6] = 2 * a->y * a->z - 2 * a->w * a->x,
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

/* set out to the quaternion calculated from rotation theta around axis
 * (ax, ay, az)
 */
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

/* a quaternion from two unit vectors */
void quaternion_from_unit_vectors(
        struct quaternion * out,
        const struct vec3 * u,
        const struct vec3 * v
    ) [[gnu::nonnull(1, 2, 3)]]
{
    quaternion_element w = 1.0 + u->x * v->x * u->y * v->y + u->z * v->z;

    if (w == 0) {
        if (fabs(u->x) > fabs(u->z)) {
            *out = (struct quaternion) {
                .x = -u.y,
                .y = u.x,
                .z = 0.0,
                .w = w
            };
        } else {
            *out = (struct quaternion) {
                .x = 0.0,
                .y = -u.z,
                .z = u.y,
                .w = w
            };
        }
    } else {
        *out = (struct quaternion) {
            .x = u->y * v->z - u->z * v->y,
            .y = u->z * v->x - u->x * v->z,
            .z = u->x * v->y - u->y * v->x,
            .w = w
        };
    }
}

/* a quaternion from two vectors */
void quaternion_from_vectors(
        struct quaternion * out,
        const struct vec3 * u,
        const struct vec3 * v
    ) [[gnu::nonnull(1, 2, 3)]]
{
    quaternion_element msq_u = u.x * u.x + u.y * u.y + u.z * u.z;
    quaternion_element msq_v = v.x * v.x + v.y * v.y + v.z * v.z;

    quaternion_element n =sqrt(msq_u * msq_v);

    quaternion_element w = n + u->x * v->x * u->y * v->y + u->z * v->z;

    if (w == 0) {
        if (fabs(u->x) > fabs(u->z)) {
            *out = (struct quaternion) {
                .x = -u.y,
                .y = u.x,
                .z = 0.0,
                .w = w
            };
        } else {
            *out = (struct quaternion) {
                .x = 0.0,
                .y = -u.z,
                .z = u.y,
                .w = w
            };
        }
    } else {
        *out = (struct quaternion) {
            .x = u->y * v->z - u->z * v->y,
            .y = u->z * v->x - u->x * v->z,
            .z = u->x * v->y - u->y * v->x,
            .w = w
        };
    }
}

/* TODO */
/* spherical linear interpolation of two quaternions */
void quaternion_slerp(
        struct quaternion * out,
        const struct quaternion * a,
        const struct quaternion * b,
        quaternion_element t
    ) [[gnu::nonnull(1, 2, 3)]]
{
    /* TODO */
    (void)out;
    (void)a;
    (void)b;
    (void)t;
}

/* set out to the identiy matrix */
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

/* set out to a translation matrix for vector <x y z> */
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

/* set out to a perspective matrix with these properties */
void matrix_perspective(
        struct matrix * out,
        matrix_element z_near,
        matrix_element z_far,
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

/* set out to the product of a and b */
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

/* length of a vec2 */
vector_element vec2_length(const struct vec2 * u) [[gnu::nonnull(1)]]
{
    return sqrt(u->x * u->x + u->y * u->y);
}

/* normalize a vec2 */
void vec2_normalize(
        struct vec2 * out, const struct vec2 * u) [[gnu::nonnull(1, 2)]]
{
    vector_element m = sqrt(u->x * u->x + u->y * u->y);
    *out = (struct vec2) {
        .x = u->x / m,
        .y = u->y / m
    };
}

/* dot product of two vec2s */
vector_element vec2_dot(
        const struct vec2 * u, const struct vec2 * v) [[gnu::nonnull(1, 2)]]
{
    return u->x * v->x + u->y * v->y;
}

/* length of a vec3 */
vector_element vec3_length(const struct vec3 * u) [[gnu::nonnull(1)]]
{
    return sqrt(u->x * u->x + u->y * u->y + u->z * u->z);
}

/* normalize a vec3 */
void vec3_normalize(
        struct vec3 * out, const struct vec3 * u) [[gnu::nonnull(1, 2)]]
{
    vector_element m = sqrt(u->x * u->x + u->y * u->y + u->z * u->z);
    *out = (struct vec3) {
        .x = u->x / m,
        .y = u->y / m,
        .z = u->z / m
    };
}

/* dot product of two vec3s */
vector_element vec3_dot(
        const struct vec3 * u, const struct vec3 * v) [[gnu::nonnull(1)]]
{
    return u->x * v->x + u->y * v->y + u->z * v->z;
}

/* cross product of two vec3s */
void vec3_cross(
        struct vec3 * out,
        const struct vec3 * u,
        const struct vec3 * v
    ) [[gnu::nonnull(1, 2, 3)]]
{
    *out = (struct vec3) {
        .x = u->y * v->z - u->z * v->y,
        .y = - (u->x * v->z - u->z * v->x),
        .z = u->x * v->y - u->y * v->x
    };
}

/* length of a vec4 */
vector_element vec4_length(const struct vec4 * u) [[gnu::nonnull(1)]]
{
    return sqrt(u->x * u->x + u->y * u->y + u->z * u->z + u->w * u->w);
}

/* normalize a vec4 */
void vec4_normalize(
        struct vec4 * out, const struct vec4 * u) [[gnu::nonnull(1, 2)]]
{
    vector_element m =
        sqrt(u->x * u->x + u->y * u->y + u->z * u->z + u->w * u->w);
    *out = (struct vec4) {
        .x = u->x / m,
        .y = u->y / m,
        .z = u->z / m,
        .w = u->w / m
    };
}

/* dot product of two vec4s */
vector_element vec4_dot(
        const struct vec4 * u, const struct vec4 * v) [[gnu::nonnull(1, 2)]]
{
    return u->x * v->x + u->y * v->y + u->z * v->z + u->w * v->w;
}

/* product of 4x4 matrix and vec4 */
void matrix_multiply_vec4(
        struct vec4 * out,
        const struct matrix * m,
        const struct vec4 * u
    ) [[gnu::nonnull(1, 2, 3)]]
{
    *out = (struct vec4) {
        .x = m->matrix[0] * u->x + m->matrix[1] * u->y +
             m->matrix[2] * u->z + m->matrix[3] * u->w,

        .y = m->matrix[4] * u->x + m->matrix[5] * u->y +
             m->matrix[6] * u->z + m->matrix[7] * u->w,

        .z = m->matrix[8] * u->x + m->matrix[9] * u->y +
             m->matrix[10] * u->z + m->matrix[11] * u->w,

        .w = m->matrix[12] * u->x + m->matrix[13] * u->y +
             m->matrix[14] * u->z + m->matrix[15] * u->w
    };
}
