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

#ifndef GRAND_TOUR_H_
#define GRAND_TOUR_H_

#include "grand-tour/topography.h"

int gt_initialise(struct gt_topography * topography, double latitude,
    double longitude, const char * path, int stack_size,
    turtle_handler_cb * handler, struct turtle_datum * datum);

/* Convert a cartesian position in local frame to an ECEF one */
void gt_to_ecef(const struct gt_topography * topography, const double * local,
    int is_a_vector, double * ecef);

/* Convert a cartesian position in local frame to a geodetic one */
int gt_to_lla(const struct gt_topography * topography, const double * local,
    double * lla);

#ifndef _GT_NO_UTM
/* Convert a cartesian position in local frame to UTM coordinates */
int gt_to_utm(const struct gt_topography * topography, const double * local,
    double * utm);
#endif

/* Convert a local direction to angular coordinates */
int gt_to_angular(const struct gt_topography * topography, const double * position,
    const double * direction, double * angular);

/* Convert a cartesian position in ECEF frame to a local one */
void gt_from_ecef(const struct gt_topography * topography, const double * ecef,
    int is_a_vector, double * local);

/* Convert a geodetic position to a cartesian one in local frame */
int gt_from_lla(const struct gt_topography * topography, const double * lla,
    double * local);

#ifndef _GT_NO_UTM
/* Convert UTM coordinates to a cartesian position in local frame */
int gt_from_utm(const struct gt_topography * topography, const double * utm,
    double * local);
#endif

/* Convert angular coordinates to a local direction */
int gt_from_angular(const struct gt_topography * topography,
    const double * position, const double * angular, double * direction);

/* Get the ground altitude in local frame coordinates or in geodetic ones */
int gt_ground_altitude(const struct gt_topography * topography,
    const double * position, int geodetic, double * altitude);

/* Compute the distance to the topography for the given point, or segment */
int gt_ground_distance(struct gt_topography * topography, const double * position,
    const double * direction, double limit, double * distance);

/* Get the normal to the ground in local or geodetic coordinates */
int gt_ground_normal(const struct gt_topography * topography,
    const double * position, int geodetic, double * normal, double * angles);

/* Check if the given local position is above the ground */
int gt_ground_above(
    const struct gt_topography * topography, const double * position, int * status);

#endif
