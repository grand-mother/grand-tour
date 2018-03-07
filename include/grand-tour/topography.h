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

#ifndef GT_TOPOGRAPHY_H_
#define GT_TOPOGRAPHY_H_

#include "turtle.h"

/* Structure for managing the topography data */
struct gt_topography {
        /* Frame origin. */
        double latitude;
        double longitude;

        /* TURTLE handles. */
        struct turtle_datum * datum;
        struct turtle_projection * projection;

        /* Flag for a flat topography. */
        int flat;
        double flat_size;

        /* Local frame parameters. */
        double origin[3];
        double base[3][3];

        /* Flag to check for TURTLE's error handling. */
        int catch;

        /* TURTLE's error handler */
        turtle_handler_cb * handler;
};

#endif
