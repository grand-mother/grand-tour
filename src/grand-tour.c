/*
 * Copyright (C) 2017 Université Clermont Auvergne, CNRS/IN2P3, LPC
 * Author: Valentin NIESS (niess@in2p3.fr)
 *
 * This software is a C99 library providing topographic utilities for GRAND. It
 * relies on the TURTLE library.
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>
 */

#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#include "grand-tour.h"

int gt_initialise(struct gt_topography * topography, double latitude,
    double longitude, const char * path, int stack_size,
    turtle_handler_cb * handler, struct turtle_datum * datum)
{
        topography->latitude = latitude;
        topography->longitude = longitude;

        /* Copy the error handler. */
        topography->catch = 0;
        topography->handler = handler;

        /* Create the datum. */
        enum turtle_return rc;
        if (datum != NULL) {
                topography->datum = datum;
        } else {
                if ((rc = turtle_datum_create(path, stack_size, NULL, NULL,
                         &topography->datum)) != TURTLE_RETURN_SUCCESS)
                        return rc;
        }

        /* Check for a flat topography and parse its size whenever provided. */
        if (strncmp(path, "flat", 4) == 0) {
                topography->flat = 1;
                if (strlen(path) > 5) {
                        char * endptr;
                        topography->flat_size = 0.5 * strtod(path + 5, &endptr);
                        if (*endptr != 0) { return -1; }
                } else {
                        topography->flat_size = 1.;
                }
        } else {
                topography->flat = 0;
                topography->flat_size = 0.;
        }

#ifndef _GT_NO_UTM
        /* Create the UTM projection. */
        static char name[8];
        const int zone = ((int)((longitude + 180) / 6) % 60) + 1;
        const char hemisphere = (latitude >= 0.) ? 'N' : 'S';
        sprintf(name, "UTM %02d%c", zone, hemisphere);
        if ((rc = turtle_projection_create(name, &topography->projection)) !=
            TURTLE_RETURN_SUCCESS)
                return rc;
#else
        topography->projection = NULL;
#endif

        /*
         * Compute the local frame using GRAND conventions, i.e. x-axis goes
         * from South to North, y-axis goes from East to West.
         */
        if ((rc = turtle_datum_ecef(topography->datum, latitude, longitude, 0.,
                 topography->origin)) != TURTLE_RETURN_SUCCESS)
                return rc;
        if ((rc = turtle_datum_direction(topography->datum, latitude, longitude,
                 0., 0., &topography->base[0][0])) != TURTLE_RETURN_SUCCESS)
                return rc;
        if ((rc = turtle_datum_direction(topography->datum, latitude, longitude,
                 270., 0., &topography->base[1][0])) != TURTLE_RETURN_SUCCESS)
                return rc;
        if ((rc = turtle_datum_direction(topography->datum, latitude, longitude,
                 0., 90., &topography->base[2][0])) != TURTLE_RETURN_SUCCESS)
                return rc;

        return TURTLE_RETURN_SUCCESS;
}

/* Convert a local position / direction to an ECEF one. */
void gt_to_ecef(const struct gt_topography * topography, const double * local,
    int is_a_vector, double * ecef)
{
        if (is_a_vector) {
                ecef[0] = 0.;
                ecef[1] = 0.;
                ecef[2] = 0.;
        } else {
                ecef[0] = topography->origin[0];
                ecef[1] = topography->origin[1];
                ecef[2] = topography->origin[2];
        }
        int i, j;
        for (i = 0; i < 3; i++)
                for (j = 0; j < 3; j++)
                        ecef[i] += topography->base[j][i] * local[j];
}

/* Convert a cartesian position in local frame to a geodetic one. */
int gt_to_lla(
    const struct gt_topography * topography, const double * local, double * lla)
{
        double ecef[3];
        gt_to_ecef(topography, local, 0, ecef);
        return turtle_datum_geodetic(
            topography->datum, ecef, lla, lla + 1, lla + 2);

        return TURTLE_RETURN_SUCCESS;
}

#ifndef _GT_NO_UTM
/* Convert a cartesian position in local frame to UTM coordinates */
int gt_to_utm(
    const struct gt_topography * topography, const double * local, double * utm)
{
        double lla[3] = { 0., 0., 0. };
        int rc;
        if ((rc = gt_to_lla(topography, local, lla)) != TURTLE_RETURN_SUCCESS)
                return rc;
        utm[2] = lla[2];

        return turtle_projection_project(
            topography->projection, lla[0], lla[1], utm, utm + 1);
}
#endif

/* Convert a local direction to angular coordinates */
int gt_to_angular(const struct gt_topography * topography,
    const double * position, const double * direction, double * angular)
{
        /* Compute the geodetic coordinates. */
        int rc;
        double lla[3] = { 0., 0., 0. };
        if ((rc = gt_to_lla(topography, position, lla)) !=
            TURTLE_RETURN_SUCCESS)
                return rc;

        /* Compute the horizontal angular coordinates. */
        double ecef[3];
        gt_to_ecef(topography, direction, 1, ecef);
        double azimuth = 0., elevation = 0.;
        if ((rc = turtle_datum_horizontal(topography->datum, lla[0], lla[1],
                 ecef, &azimuth, &elevation)) != TURTLE_RETURN_SUCCESS)
                return rc;

        /* Convert to GRAND conventions and return. */
        angular[0] = 90. - elevation;
        angular[1] = -azimuth;

        return TURTLE_RETURN_SUCCESS;
}

/* Convert an ECEF position/direction to a local one. */
void gt_from_ecef(const struct gt_topography * topography, const double * ecef,
    int is_a_vector, double * local)
{
        const double eps = is_a_vector ? 0. : 1.;
        int i, j;
        for (i = 0; i < 3; i++) {
                local[i] = 0.;
                for (j = 0; j < 3; j++)
                        local[i] += topography->base[i][j] *
                            (ecef[j] - eps * topography->origin[j]);
        }
}

/* Convert a geodetic position to a cartesian one in local frame */
int gt_from_lla(
    const struct gt_topography * topography, const double * lla, double * local)
{
        enum turtle_return rc;
        double ecef[3];
        if ((rc = turtle_datum_ecef(topography->datum, lla[0], lla[1], lla[2],
                 ecef)) != TURTLE_RETURN_SUCCESS)
                return rc;
        gt_from_ecef(topography, ecef, 0, local);
        return TURTLE_RETURN_SUCCESS;
}

#ifndef _GT_NO_UTM
/* Convert UTM coordinates to a cartesian position in local frame */
int gt_from_utm(
    const struct gt_topography * topography, const double * utm, double * local)
{
        double lla[3] = { 0., 0., utm[2] };
        int rc;
        if ((rc = turtle_projection_unproject(topography->projection, utm[0],
                 utm[1], lla, lla + 1)) != TURTLE_RETURN_SUCCESS)
                return rc;
        return gt_from_lla(topography, lla, local);
}
#endif

/* Convert angular coordinates to a local direction */
int gt_from_angular(const struct gt_topography * topography,
    const double * position, const double * angular, double * direction)
{
        /* Compute the geodetic coordinates. */
        int rc;
        double lla[3] = { 0., 0., 0. };
        if ((rc = gt_to_lla(topography, position, lla)) !=
            TURTLE_RETURN_SUCCESS)
                return rc;

        /* Compute the local direction. */
        double ecef[3];
        const double azimuth = -angular[1];
        const double elevation = 90. - angular[0];
        if ((rc = turtle_datum_direction(topography->datum, lla[0], lla[1],
                 azimuth, elevation, ecef)) != TURTLE_RETURN_SUCCESS)
                return rc;
        gt_from_ecef(topography, ecef, 1, direction);

        return TURTLE_RETURN_SUCCESS;
}

/* Encapsulate TURTLE calls to elevation in order to check for a flat
 * topography.
 */
static enum turtle_return ground_elevation(
    const struct gt_topography * topography, double latitude, double longitude,
    double * z)
{
        if (topography->flat) {
                *z = 0.;
                if ((fabs(latitude - topography->latitude) >
                        topography->flat_size) ||
                    (fabs(longitude - topography->longitude) >
                        topography->flat_size)) {
                        if (!topography->catch && (topography->handler != NULL))
                                topography->handler(TURTLE_RETURN_PATH_ERROR,
                                    (turtle_caller_t *)turtle_datum_elevation);
                        return TURTLE_RETURN_PATH_ERROR;
                } else {
                        return TURTLE_RETURN_SUCCESS;
                }
        } else {
                return turtle_datum_elevation(
                    topography->datum, latitude, longitude, z);
        }
}

/* Estimate the ground elevation in the local frame. */
static enum turtle_return ground_elevation_local(
    const struct gt_topography * topography, double x, double y, double * zg)
{
        /* Let us compute the altitude in local coordinates using an iterative
         * method. We need to find the ground position that projects at (x, y).
         * Let us Initialise the search on the ellipsoid, i.e. at zero altitude.
         */
        double local[3] = { x, y, 0. };
        double new[3] = { 0., 0., 0. };
        int _;
        for (_ = 0; _ < 5; _++) {
                /* Get the ground altitude for the current guess. */
                enum turtle_return rc;
                double lla[3] = { 0., 0., 0. };
                if ((rc = gt_to_lla(topography, local, lla)) !=
                    TURTLE_RETURN_SUCCESS)
                        return rc;
                if ((rc = ground_elevation(topography, lla[0], lla[1],
                         lla + 2)) != TURTLE_RETURN_SUCCESS)
                        return rc;

                /* Compute the corresponding local coordinates. */
                if ((rc = gt_from_lla(topography, lla, new)) !=
                    TURTLE_RETURN_SUCCESS)
                        return rc;

                /* Check for convergence and update. */
                const double dx = x - new[0];
                const double dy = y - new[1];
                if ((fabs(dx) < 1E-03) && (fabs(dy) < 1E-03)) break;
                local[0] += dx;
                local[1] += dy;
        }
        *zg = new[2];

        return TURTLE_RETURN_SUCCESS;
}

/* Get the ground altitude in local frame coordinates or in geodetic ones */
int gt_ground_altitude(const struct gt_topography * topography,
    const double * position, int geodetic, double * altitude)
{
        if (geodetic) {
                return ground_elevation(
                    topography, position[0], position[1], altitude);
        } else {
                return ground_elevation_local(
                    topography, position[0], position[1], altitude);
        }
}

/* Compute the vertical distance to the ground. */
static enum turtle_return vertical_distance(
    const struct gt_topography * topography, const double * local,
    double * distance)
{
        enum turtle_return rc;
        double lla[3];
        if ((rc = gt_to_lla(topography, local, lla)) != TURTLE_RETURN_SUCCESS)
                return rc;
        double z;
        if ((rc = ground_elevation(topography, lla[0], lla[1], &z)) !=
            TURTLE_RETURN_SUCCESS)
                return rc;

        *distance = lla[2] - z;
        return TURTLE_RETURN_SUCCESS;
}

/* Check if the given local position is above the ground */
int gt_ground_above(const struct gt_topography * topography,
    const double * position, int * status)
{
        int rc;
        double distance;
        if ((rc = vertical_distance(topography, position, &distance)) !=
            TURTLE_RETURN_SUCCESS)
                return rc;
        *status = (distance > 0.);
        return TURTLE_RETURN_SUCCESS;
}

/* Compute the distance to the topography for the given point, or segment */
int gt_ground_distance(struct gt_topography * topography,
    const double * position, const double * direction, double limit,
    double * distance_ptr)
{
        if (direction == NULL) {
                /* This is a point. Let us return the vertical distance. */
                return vertical_distance(topography, position, distance_ptr);
        }

        /* Let us step along the given direction until the ground is
         * reached.
         */
        const double resolution = 1E-02;
        int rc, above;
        if ((rc = gt_ground_above(topography, position, &above)) !=
            TURTLE_RETURN_SUCCESS)
                return rc;
        double r0[3] = { position[0], position[1], position[2] };
        double distance = 0.;
        turtle_handler(NULL);
        topography->catch = 1;
        for (;;) {
                double step;
                if (vertical_distance(topography, r0, &step) !=
                    TURTLE_RETURN_SUCCESS)
                        break;
                step = 0.5 * fabs(step);
                if (step < resolution) step = resolution;
                double r1[3] = { r0[0] + step * direction[0],
                        r0[1] + step * direction[1],
                        r0[2] + step * direction[2] };
                int above1;
                if (gt_ground_above(topography, r1, &above1) !=
                    TURTLE_RETURN_SUCCESS)
                        break;
                if (above1 != above) {
                        /* A medium change was found. Let us refine it with a
                         * binary search.
                         */
                        while (step > resolution) {
                                double r2[3] = { 0.5 * (r0[0] + r1[0]),
                                        0.5 * (r0[1] + r1[1]),
                                        0.5 * (r0[2] + r1[2]) };
                                int above2 = 0;
                                gt_ground_above(topography, r2, &above2);
                                if (above2 == above)
                                        memcpy(r0, r2, sizeof(r0));
                                else
                                        memcpy(r1, r2, sizeof(r1));
                                step *= 0.5;
                        }
                        distance = 0.;
                        int i;
                        for (i = 0; i < 3; i++) {
                                const double d = r1[i] - position[i];
                                distance += d * d;
                        }
                        distance = (distance > 0.) ? sqrt(distance) : 0.;
                        if ((limit > 0.) && (distance > limit)) break;
                        if (!above) distance = -distance;
                        turtle_handler(topography->handler);
                        topography->catch = 0;
                        *distance_ptr = distance;
                        return TURTLE_RETURN_SUCCESS;
                } else if (fabs(r1[2] > 1E+04))
                        break;
                if (limit > 0.) {
                        distance += step;
                        if (distance > limit) break;
                }
                memcpy(r0, r1, sizeof(r0));
        }

        turtle_handler(topography->handler);
        topography->catch = 0;
        *distance_ptr = -1.;
        return TURTLE_RETURN_SUCCESS;
}

/* Get the normal to the ground in local frame coordinates */
static int ground_normal_local(const struct gt_topography * topography,
    const double * position, double step, double * normal)
{
        /* Estimate the local slope using an adaptation of S. Le Coz's
         * algorithm.
         */
        enum turtle_return rc;
        const double ds = (step <= 0.) ? 10. : step;
        double zx0, zx1, zy0, zy1;
        if ((rc = ground_elevation_local(topography, position[0] - 0.5 * ds,
                 position[1], &zx0)) != TURTLE_RETURN_SUCCESS)
                return rc;
        if ((rc = ground_elevation_local(topography, position[0] + 0.5 * ds,
                 position[1], &zx1)) != TURTLE_RETURN_SUCCESS)
                return rc;
        if ((rc = ground_elevation_local(topography, position[0],
                 position[1] - 0.5 * ds, &zy0)) != TURTLE_RETURN_SUCCESS)
                return rc;
        if ((rc = ground_elevation_local(topography, position[0],
                 position[1] + 0.5 * ds, &zy1)) != TURTLE_RETURN_SUCCESS)
                return rc;

        double u[3] = { ds, 0., zx1 - zx0 };
        double v[3] = { 0., ds, zy1 - zy0 };
        normal[0] = u[1] * v[2] - u[2] * v[1];
        normal[1] = u[2] * v[0] - u[0] * v[2];
        normal[2] = u[0] * v[1] - u[1] * v[0];
        double d = normal[0] * normal[0] + normal[1] * normal[1] +
            normal[2] * normal[2];
        if (d <= FLT_EPSILON)
                d = 0.;
        else
                d = 1. / sqrt(d);
        normal[0] *= d;
        normal[1] *= d;
        normal[2] *= d;

        return TURTLE_RETURN_SUCCESS;
}

/* Get the normal to the ground in local or geodetic coordinates */
int gt_ground_normal(const struct gt_topography * topography,
    const double * position, int geodetic, double step, double * normal,
    double * angles)
{
        int rc;
        double tmp[3];
        if (topography->flat) {
                /* Compute the local verticale. */
                if (!geodetic) {
                        double z;
                        if ((rc = ground_elevation_local(
                                 topography, position[0], position[1], &z)) !=
                            TURTLE_RETURN_SUCCESS)
                                return rc;

                        double local[3] = { position[0], position[1], z };
                        gt_to_lla(topography, local, tmp);
                        position = tmp;
                }

                double ecef[3];
                if ((rc = turtle_datum_direction(topography->datum, position[0],
                         position[1], 0., 90., ecef)) != TURTLE_RETURN_SUCCESS)
                        return rc;
                gt_from_ecef(topography, ecef, 1, normal);
        } else {
                /* Estimate the local slope of the ground. */
                if (geodetic) {
                        double z;
                        if ((rc = ground_elevation(topography, position[0],
                                 position[1], &z)) != TURTLE_RETURN_SUCCESS)
                                return rc;

                        double lla[3] = { position[0], position[1], z };
                        if ((rc = gt_from_lla(topography, lla, tmp)) !=
                            TURTLE_RETURN_SUCCESS)
                                return rc;
                        position = tmp;
                }

                if ((rc = ground_normal_local(topography, position, step,
                         normal)) != TURTLE_RETURN_SUCCESS)
                        return rc;
        }

        if (angles != NULL) {
                angles[1] = atan2(normal[1], normal[0]) * 180. / M_PI;
                angles[0] = acos(normal[2]) * 180. / M_PI;
        }

        return TURTLE_RETURN_SUCCESS;
}
