/* File: src/test.c
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
#include <stdio.h>
#include <stdlib.h>

size_t unexpected = 0;
double accuracy = 1e-12;

char ** string_pool = NULL;
size_t string_pool_size = 0;

void free_string_pool()
{
    for (size_t i = 0; i < string_pool_size; i++) {
        free(string_pool[i]);
    }
    free(string_pool);
    string_pool = NULL;
    string_pool_size = 0;
}

const char * strof_q(const struct quaternion * q)
{
    char * buffer = malloc(1024);
    snprintf(buffer, 1024, "<%f %f %f %f>", q->x, q->y, q->z, q->w);
    string_pool = realloc(
            string_pool, sizeof(*string_pool) * (string_pool_size + 1));
    string_pool[string_pool_size++] = buffer;
    return buffer;
}

const char * strof_v2(const struct vec2 * u)
{
    char * buffer = malloc(1024);
    snprintf(buffer, 1024, "<%f %f>", u->x, u->y);
    string_pool = realloc(
            string_pool, sizeof(*string_pool) * (string_pool_size + 1));
    string_pool[string_pool_size++] = buffer;
    return buffer;
}

const char * strof_v3(const struct vec3 * u)
{
    char * buffer = malloc(1024);
    snprintf(buffer, 1024, "<%f %f %f>", u->x, u->y, u->z);
    string_pool = realloc(
            string_pool, sizeof(*string_pool) * (string_pool_size + 1));
    string_pool[string_pool_size++] = buffer;
    return buffer;
}

const char * strof_v4(const struct vec4 * u)
{
    char * buffer = malloc(1024);
    snprintf(buffer, 1024, "<%f %f %f %f>", u->x, u->y, u->z, u->w);
    string_pool = realloc(
            string_pool, sizeof(*string_pool) * (string_pool_size + 1));
    string_pool[string_pool_size++] = buffer;
    return buffer;
}

const char * strof_m(const struct matrix * m)
{
    char * buffer = malloc(1024);
    snprintf(
            buffer,
            1024,
            "[[%f %f %f %f]\n"
            " [%f %f %f %f]\n"
            " [%f %f %f %f]\n"
            " [%f %f %f %f]]",
            m->matrix[0],
            m->matrix[1],
            m->matrix[2],
            m->matrix[3],
            m->matrix[4],
            m->matrix[5],
            m->matrix[6],
            m->matrix[7],
            m->matrix[8],
            m->matrix[9],
            m->matrix[10],
            m->matrix[11],
            m->matrix[12],
            m->matrix[13],
            m->matrix[14],
            m->matrix[15]
        );
    string_pool = realloc(
            string_pool, sizeof(*string_pool) * (string_pool_size + 1));
    string_pool[string_pool_size++] = buffer;
    return buffer;
}

void expect_vec3(const struct vec3 * a, const struct vec3 * b)
{
    if (fabs(a->x - b->x) > accuracy ||
        fabs(a->y - b->y) > accuracy ||
        fabs(a->z - b->z) > accuracy) {
        printf("UNEXPECTED: %s != %s\n", strof_v3(a), strof_v3(b));
        unexpected++;
    }
}

void expect_q(const struct quaternion * a, const struct quaternion * b)
{
    if (fabs(a->x - b->x) > accuracy ||
        fabs(a->y - b->y) > accuracy ||
        fabs(a->z - b->z) > accuracy ||
        fabs(a->w - b->w) > accuracy) {
        printf("UNEXPECTED: %s != %s\n", strof_q(a), strof_q(b));
        unexpected++;
    }
}

static void perspective_test(
        const char * name,
        struct vec4 * in,
        struct matrix * matrix
    );

int main(int argc, char ** argv)
{
    (void)argc;
    (void)argv;

    struct quaternion a;
    struct matrix m;

    {
        struct vec4 u = { .x = 0.0, .y = 0.0, .z = 1.0, .w = 1.0 };
        printf("input = %s\n", strof_v4(&u));
        quaternion_from_axis_angle(&a, 1.0, 0.0, 0.0, M_PI / 2);
        printf("rotate by x axis, pi / 2 (90°)\n");
        printf("quaternion = %s\n", strof_q(&a));
        quaternion_matrix(&m, &a);
        printf("as matrix =\n%s\n", strof_m(&m));
        matrix_multiply_vec4(&u, &m, &u);
        printf("as vec4 = %s\n", strof_v4(&u));
        struct vec3 v = { .x = u.x / u.w, .y = u.y / u.w, .z = u.z / u.w };
        printf("output = %s\n", strof_v3(&v));
        expect_vec3(&v, &(struct vec3){0, -1, 0});
    }
    printf("\n");

    {
        struct vec4 u = { .x = 0.0, .y = 0.0, .z = 1.0, .w = 1.0 };
        printf("input = %s\n", strof_v4(&u));
        quaternion_from_axis_angle(&a, 1.0, 0.0, 0.0, -M_PI / 2);
        printf("rotate by x axis, -pi / 2 (-90°)\n");
        printf("quaternion = %s\n", strof_q(&a));
        quaternion_matrix(&m, &a);
        printf("as matrix =\n%s\n", strof_m(&m));
        matrix_multiply_vec4(&u, &m, &u);
        printf("as vec4 = %s\n", strof_v4(&u));
        struct vec3 v = { .x = u.x / u.w, .y = u.y / u.w, .z = u.z / u.w };
        printf("output = %s\n", strof_v3(&v));
        expect_vec3(&v, &(struct vec3){0, 1, 0});
    }
    printf("\n");

    {
        struct vec4 u = { .x = 0.0, .y = 0.0, .z = 1.0, .w = 1.0 };
        printf("input = %s\n", strof_v4(&u));
        quaternion_from_axis_angle(&a, 0.0, 1.0, 0.0, M_PI / 2);
        printf("rotate by y axis, pi / 2 (90°)\n");
        printf("quaternion = %s\n", strof_q(&a));
        quaternion_matrix(&m, &a);
        printf("as matrix =\n%s\n", strof_m(&m));
        matrix_multiply_vec4(&u, &m, &u);
        printf("as vec4 = %s\n", strof_v4(&u));
        struct vec3 v = { .x = u.x / u.w, .y = u.y / u.w, .z = u.z / u.w };
        printf("output = %s\n", strof_v3(&v));
        expect_vec3(&v, &(struct vec3){1, 0, 0});
    }
    printf("\n");

    {
        struct vec4 u = { .x = 0.0, .y = 0.0, .z = 1.0, .w = 1.0 };
        printf("input = %s\n", strof_v4(&u));
        quaternion_from_axis_angle(&a, 0.0, 1.0, 0.0, -M_PI / 2);
        printf("rotate by y axis, -pi / 2 (-90°)\n");
        printf("quaternion = %s\n", strof_q(&a));
        quaternion_matrix(&m, &a);
        printf("as matrix =\n%s\n", strof_m(&m));
        matrix_multiply_vec4(&u, &m, &u);
        printf("as vec4 = %s\n", strof_v4(&u));
        struct vec3 v = { .x = u.x / u.w, .y = u.y / u.w, .z = u.z / u.w };
        printf("output = %s\n", strof_v3(&v));
        expect_vec3(&v, &(struct vec3){-1, 0, 0});
    }
    printf("\n");

    {
        struct vec4 u = { .x = 0.0, .y = 0.0, .z = 1.0, .w = 1.0 };
        printf("input = %s\n", strof_v4(&u));
        quaternion_from_axis_angle(&a, 0.0, 0.0, 1.0, -M_PI / 2);
        printf("rotate by z axis, -pi / 2 (-90°)\n");
        printf("quaternion = %s\n", strof_q(&a));
        quaternion_matrix(&m, &a);
        printf("as matrix =\n%s\n", strof_m(&m));
        matrix_multiply_vec4(&u, &m, &u);
        printf("as vec4 = %s\n", strof_v4(&u));
        struct vec3 v = { .x = u.x / u.w, .y = u.y / u.w, .z = u.z / u.w };
        printf("output = %s\n", strof_v3(&v));
        expect_vec3(&v, &(struct vec3){0, 0, 1});
    }
    printf("\n");

    {
        struct vec4 u = { .x = 0.0, .y = 0.0, .z = 1.0, .w = 1.0 };
        printf("input = %s\n", strof_v4(&u));
        quaternion_from_axis_angle(&a, 0.0, 0.0, 1.0, M_PI / 2);
        printf("rotate by z axis, pi / 2 (90°)\n");
        printf("quaternion = %s\n", strof_q(&a));
        quaternion_matrix(&m, &a);
        printf("as matrix =\n%s\n", strof_m(&m));
        matrix_multiply_vec4(&u, &m, &u);
        printf("as vec4 = %s\n", strof_v4(&u));
        struct vec3 v = { .x = u.x / u.w, .y = u.y / u.w, .z = u.z / u.w };
        printf("output = %s\n", strof_v3(&v));
        expect_vec3(&v, &(struct vec3){0, 0, 1});
    }
    printf("\n");

    {
        struct vec4 u = { .x = 0.0, .y = 0.0, .z = 1.0, .w = 1.0 };
        printf("input = %s\n", strof_v4(&u));
        quaternion_from_axis_angle(&a, 1.0, 0.0, 0.0, M_PI);
        printf("rotate by x axis, pi (180°)\n");
        printf("quaternion = %s\n", strof_q(&a));
        quaternion_matrix(&m, &a);
        printf("as matrix =\n%s\n", strof_m(&m));
        matrix_multiply_vec4(&u, &m, &u);
        printf("as vec4 = %s\n", strof_v4(&u));
        struct vec3 v = { .x = u.x / u.w, .y = u.y / u.w, .z = u.z / u.w };
        printf("output = %s\n", strof_v3(&v));
        expect_vec3(&v, &(struct vec3){0, 0, -1});
    }
    printf("\n");

    {
        struct vec4 u = { .x = 0.0, .y = 0.0, .z = 1.0, .w = 1.0 };
        printf("input = %s\n", strof_v4(&u));
        quaternion_from_axis_angle(&a, 0.0, 1.0, 0.0, M_PI);
        printf("rotate by y axis, pi (180°)\n");
        printf("quaternion = %s\n", strof_q(&a));
        quaternion_matrix(&m, &a);
        printf("as matrix =\n%s\n", strof_m(&m));
        matrix_multiply_vec4(&u, &m, &u);
        printf("as vec4 = %s\n", strof_v4(&u));
        struct vec3 v = { .x = u.x / u.w, .y = u.y / u.w, .z = u.z / u.w };
        printf("output = %s\n", strof_v3(&v));
        expect_vec3(&v, &(struct vec3){0, 0, -1});
    }
    printf("\n");

    {
        struct vec4 u = { .x = 0.0, .y = 0.0, .z = 1.0, .w = 1.0 };
        printf("input = %s\n", strof_v4(&u));
        quaternion_from_axis_angle(&a, 0.0, 1.0, 0.0, M_PI / 10);
        printf("rotate by x axis, pi / 10 x 10\n");
        printf("quaternion a = %s\n", strof_q(&a));
        struct quaternion b;
        quaternion_identity(&b);
        for (size_t i = 0; i < 10; i++) {
            quaternion_multiply(&b, &b, &a);
        }
        quaternion_normalize(&b, &b);
        printf("quaternion b = %s\n", strof_q(&b));
        quaternion_matrix(&m, &b);
        printf("as matrix =\n%s\n", strof_m(&m));
        matrix_multiply_vec4(&u, &m, &u);
        printf("as vec4 = %s\n", strof_v4(&u));
        struct vec3 v = { .x = u.x / u.w, .y = u.y / u.w, .z = u.z / u.w };
        printf("output = %s\n", strof_v3(&v));
        expect_vec3(&v, &(struct vec3){0, 0, -1});
    }
    printf("\n");

    {
        struct vec4 u = { .x = 0.0, .y = 0.0, .z = 1.0, .w = 1.0 };
        printf("input = %s\n", strof_v4(&u));
        quaternion_from_axis_angle(&a, 1.0, 0.0, 0.0, M_PI / 10);
        printf("rotate by x axis, pi / 10 x 10\n");
        printf("quaternion a = %s\n", strof_q(&a));
        struct quaternion b;
        quaternion_identity(&b);
        for (size_t i = 0; i < 10; i++) {
            quaternion_multiply(&b, &b, &a);
        }
        quaternion_normalize(&b, &b);
        printf("quaternion b = %s\n", strof_q(&b));
        quaternion_matrix(&m, &b);
        printf("as matrix =\n%s\n", strof_m(&m));
        matrix_multiply_vec4(&u, &m, &u);
        printf("as vec4 = %s\n", strof_v4(&u));
        struct vec3 v = { .x = u.x / u.w, .y = u.y / u.w, .z = u.z / u.w };
        printf("output = %s\n", strof_v3(&v));
        expect_vec3(&v, &(struct vec3){0, 0, -1});
    }
    printf("\n");

    {
        struct quaternion a = { 1.0, 2.0, 4.0, 3.0 };
        struct quaternion b = { 3.0, 2.0, 3.0, 2.0 };
        struct quaternion out;
        quaternion_multiply(&out, &a, &b);
        printf("output = %s\n", strof_q(&out));
        expect_q(&out, &(struct quaternion) { 9.0, 19.0, 13.0, -13.0 });
    }
    printf("\n");

    {
        struct quaternion a = { 2.0, 1.0, 7.0, 4.0 };
        struct quaternion b = { 2.0, -3.0, 4.0, 11.0 };
        struct quaternion out;
        quaternion_multiply(&out, &a, &b);
        printf("output = %s\n", strof_q(&out));
        expect_q(&out, &(struct quaternion) { 55.0, 5.0, 85.0, 15.0 });
    }
    printf("\n");

    {
        /* in world coordinates:
         *   -x is left, x is right
         *   -y is down, y is up
         *   -z is farther, z is nearer
         */

        struct vec4 left = { -1.0, 0.0, 0.0, 1.0 };
        struct vec4 right = { 1.0, 0.0, 0.0, 1.0 };
        struct vec4 down = { 0.0, -1.0, 0.0, 1.0 };
        struct vec4 up = { 0.0, 1.0, 0.0, 1.0 };
        struct vec4 far = { 0.0, 0.0, -1.0, 1.0 };
        struct vec4 near = { 0.0, 0.0, 1.0, 1.0 };
        struct vec4 center = { 0.0, 0.0, -0.0, 1.0 };

        /* in vulkan:
         * y goes from (-1, +1) from top to bottom
         * x goes from (-1, +1) from left to right
         * z goes from (0, 1) from near to far
         *
         * we must also correct for aspect (ratio of screen x to y)
         */

        struct matrix matrix_a;
        struct matrix matrix;
        matrix_translation(&matrix_a, 0, 0, -5);
        matrix_perspective(&matrix, -0.1, -10.0, M_PI / 4, 1);
        matrix_multiply(&matrix, &matrix, &matrix_a);

        perspective_test("left", &left, &matrix);
        perspective_test("right", &right, &matrix);
        perspective_test("down", &down, &matrix);
        perspective_test("up", &up, &matrix);
        perspective_test("far", &far, &matrix);
        perspective_test("near", &near, &matrix);
        perspective_test("center", &center, &matrix);
    }
    printf("\n");


    printf("%zu unexpected\n", unexpected);

    free_string_pool();
}

static void perspective_test(
        const char * name,
        struct vec4 * in,
        struct matrix * matrix
    )
{
    struct vec4 out;
    matrix_multiply_vec4(&out, matrix, in);
    struct vec3 ndc = { out.x / out.w, out.y / out.w, out.z / out.w };
    printf("%s %s\n", name, strof_v3(&ndc));
}
