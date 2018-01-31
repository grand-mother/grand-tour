/*
 * Copyright (C) 2017 Université Clermont Auvergne, CNRS/IN2P3, LPC
 * Author: Valentin NIESS (niess@in2p3.fr)
 *
 * This software is a Python C extension module providing topographic
 * utilities for GRAND. It relies on the TURTLE library.
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

#include "Python.h"
#include "turtle.h"

/* Python exception for TURTLE errors. */
static PyObject * TurtleError = NULL;

/* Raise a Python exception whenever a TURTLE library error is triggered. */
static void handle_turtle_error(enum turtle_return rc, turtle_caller_t * caller)
{
        static char msg[256];
        sprintf(msg, "%s [%s]", turtle_strfunc(caller), turtle_strerror(rc));
        PyErr_SetString(TurtleError, msg);
}

/* The low level topography data. */
typedef struct {
        PyObject_HEAD
            /* TURTLE handles. */
            struct turtle_datum * datum;
        struct turtle_projection * projection;
        /* Flag to check for TURTLE's error handling. */
        int catch;
        /* Flag for a flat topography. */
        int flat;
        double flat_size;
        /* Local frame parameters. */
        double latitude;
        double longitude;
        double origin[3];
        double base[3][3];
} TopographyObject;

/* Initialise a Topography object. */
static int topography_init(
    TopographyObject * self, PyObject * args, PyObject * kwargs)
{
        /* Parse the arguments. */
        double latitude, longitude;
        char * path;
        int stack_size = 4;

        static char * kwlist[] = { "latitude", "longitude", "path",
                "stack_size", NULL };
        if (!PyArg_ParseTupleAndKeywords(args, kwargs, "dds|i", kwlist,
                &latitude, &longitude, &path, &stack_size))
                return -1;
        self->catch = 0;
        self->latitude = latitude;
        self->longitude = longitude;

        /* Create the datum. */
        if ((turtle_datum_create(path, stack_size, NULL, NULL, &self->datum)) !=
            TURTLE_RETURN_SUCCESS)
                return -1;

        /* Check for a flat topography and parse its size whenever provided. */
        if (strncmp(path, "flat", 4) == 0) {
                self->flat = 1;
                if (strlen(path) > 5) {
                        char * endptr;
                        self->flat_size = 0.5 * strtod(path + 5, &endptr);
                        if (*endptr != 0) {
                                PyErr_SetString(PyExc_ValueError,
                                    "invalid size for flat topography");
                                return -1;
                        }
                } else {
                        self->flat_size = 1.;
                }
        } else {
                self->flat = 0;
                self->flat_size = 0.;
        }

        /* Create the UTM projection. */
        static char name[8];
        const int zone = ((int)((longitude + 180) / 6) % 60) + 1;
        const char hemisphere = (latitude >= 0.) ? 'N' : 'S';
        sprintf(name, "UTM %02d%c", zone, hemisphere);
        if ((turtle_projection_create(name, &self->projection)) !=
            TURTLE_RETURN_SUCCESS)
                return -1;

        /*
         * Compute the local frame using GRAND conventions, i.e. x-axis goes
         * from South to North, y-axis goes from East to West.
         */
        if ((turtle_datum_ecef(self->datum, latitude, longitude, 0.,
                self->origin)) != TURTLE_RETURN_SUCCESS)
                return -1;
        if ((turtle_datum_direction(self->datum, latitude, longitude, 0., 0.,
                &self->base[0][0])) != TURTLE_RETURN_SUCCESS)
                return -1;
        if ((turtle_datum_direction(self->datum, latitude, longitude, 270., 0.,
                &self->base[1][0])) != TURTLE_RETURN_SUCCESS)
                return -1;
        if ((turtle_datum_direction(self->datum, latitude, longitude, 0., 90.,
                &self->base[2][0])) != TURTLE_RETURN_SUCCESS)
                return -1;

        return 0;
}

/* Free a topography object. */
static void topography_destroy(TopographyObject * self)
{
        turtle_datum_destroy(&self->datum);
        turtle_projection_destroy(&self->projection);
        Py_TYPE(self)->tp_free((PyObject *)self);
}

/* Utility function for parsing an arbitrary sequence as a vector. */
static int parse_vector(PyObject * object, int n, double * vector)
{
        const int size = PySequence_Length(object);
        if (size < n) return -1;
        int i;
        for (i = 0; i < n; i++) {
                PyObject * value = PySequence_GetItem(object, i);
                if (value == NULL) return -1;
                PyObject * float_ = PyNumber_Float(value);
                Py_DECREF(value);
                if (float_ == NULL) return -1;
                vector[i] = PyFloat_AsDouble(float_);
                Py_DECREF(float_);
                if (PyErr_Occurred() != NULL) return -1;
        }
        return 0;
}

/* Convert a local position/direction to an ECEF one. */
static void local_to_ecef(
    TopographyObject * self, const double * local, double * ecef, int vector)
{
        if (vector) {
                ecef[0] = 0.;
                ecef[1] = 0.;
                ecef[2] = 0.;
        } else {
                ecef[0] = self->origin[0];
                ecef[1] = self->origin[1];
                ecef[2] = self->origin[2];
        }
        int i, j;
        for (i = 0; i < 3; i++)
                for (j = 0; j < 3; j++) ecef[i] += self->base[j][i] * local[j];
}

static PyObject * topography_local_to_ecef(
        TopographyObject * self, PyObject * args, PyObject * kwargs)
{
        /* Parse the arguments. */
        PyObject * local_obj = NULL;
        int vector = 0;

        static char * kwlist[] = { "local", "vector", NULL };
        if (!PyArg_ParseTupleAndKeywords(
            args, kwargs, "O|b", kwlist, &local_obj, &vector))
            return NULL;

        double local[3];
        if (parse_vector(local_obj, 3, local) != 0) return NULL;
        double ecef[3];
        local_to_ecef(self, local, ecef, vector);

        return Py_BuildValue("(d,d,d)", ecef[0], ecef[1], ecef[2]);
}

/* Convert a cartesian position in local frame to a geodetic one. */
static int local_to_lla(TopographyObject * self, double * local,
    double * latitude, double * longitude, double * altitude)
{
        double ecef[3];
        local_to_ecef(self, local, ecef, 0);
        return turtle_datum_geodetic(
            self->datum, ecef, latitude, longitude, altitude);
}

static PyObject * topography_local_to_lla(
    TopographyObject * self, PyObject * position)
{
        double local[3];
        if (parse_vector(position, 3, local) != 0) return NULL;
        double latitude = 0, longitude = 0, altitude = 0;
        if (local_to_lla(self, local, &latitude, &longitude, &altitude) != 0)
                return NULL;
        return Py_BuildValue("ddd", latitude, longitude, altitude);
}

/* Convert a cartesian position in local frame to UTM coordinates. */
static PyObject * topography_local_to_utm(
    TopographyObject * self, PyObject * position)
{
        double local[3];
        if (parse_vector(position, 3, local) != 0) return NULL;
        double latitude = 0., longitude = 0., altitude = 0.;
        if (local_to_lla(self, local, &latitude, &longitude, &altitude) != 0)
                return NULL;
        double x = 0., y = 0.;
        if (turtle_projection_project(self->projection, latitude, longitude, &x,
                &y) != TURTLE_RETURN_SUCCESS)
                return NULL;
        return Py_BuildValue("ddd", x, y, altitude);
}

/* Convert a local direction to GRAND angular coordinates. */
static PyObject * topography_local_to_angular(
    TopographyObject * self, PyObject * args)
{
        /* Parse the arguments. */
        PyObject *position_obj = NULL, *direction_obj = NULL;
        if (!PyArg_ParseTuple(args, "OO", &position_obj, &direction_obj))
                return NULL;
        double local[3];
        if (parse_vector(position_obj, 3, local) != 0) return NULL;
        double direction[3];
        if (parse_vector(direction_obj, 3, direction) != 0) return NULL;

        /* Compute the geodetic coordinates. */
        double latitude = 0., longitude = 0., altitude = 0.;
        if (local_to_lla(self, local, &latitude, &longitude, &altitude) != 0)
                return NULL;

        /* Compute the horizontal angular coordinates. */
        double ecef[3];
        local_to_ecef(self, direction, ecef, 1);
        double azimuth = 0., elevation = 0.;
        if (turtle_datum_horizontal(self->datum, latitude, longitude, ecef,
                &azimuth, &elevation) != TURTLE_RETURN_SUCCESS)
                return NULL;

        /* Convert to GRAND conventions and return. */
        const double theta = 90. - elevation;
        const double phi = -azimuth;
        return Py_BuildValue("dd", theta, phi);
}

/* Convert an ECEF position/direction to a local one. */
static void ecef_to_local(
    TopographyObject * self, const double * ecef, double * local, int vector)
{
        const double eps = vector ? 0. : 1.;
        int i, j;
        for (i = 0; i < 3; i++) {
                local[i] = 0.;
                for (j = 0; j < 3; j++)
                        local[i] += self->base[i][j] * (ecef[j] -
                            eps * self->origin[j]);
        }
}

static PyObject * topography_ecef_to_local(
        TopographyObject * self, PyObject * args, PyObject * kwargs)
{
        /* Parse the arguments. */
        PyObject * ecef_obj = NULL;
        int vector = 0;

        static char * kwlist[] = { "ecef", "vector", NULL };
        if (!PyArg_ParseTupleAndKeywords(
            args, kwargs, "O|b", kwlist, &ecef_obj, &vector))
            return NULL;

        double ecef[3];
        if (parse_vector(ecef_obj, 3, ecef) != 0) return NULL;
        double local[3];
        ecef_to_local(self, ecef, local, vector);

        return Py_BuildValue("(d,d,d)", local[0], local[1], local[2]);
}

/* Convert a geodetic position to a cartesian one in local frame. */
static enum turtle_return lla_to_local(TopographyObject * self, double latitude,
    double longitude, double altitude, double * local)
{
        enum turtle_return rc;
        double ecef[3];
        if ((rc = turtle_datum_ecef(self->datum, latitude, longitude, altitude,
                 ecef)) != TURTLE_RETURN_SUCCESS)
                return rc;
        ecef_to_local(self, ecef, local, 0);
        return TURTLE_RETURN_SUCCESS;
}

static PyObject * topography_lla_to_local(
    TopographyObject * self, PyObject * args)
{
        double latitude, longitude, altitude;
        if (!PyArg_ParseTuple(args, "ddd", &latitude, &longitude, &altitude))
                return NULL;
        double local[3] = { 0., 0., 0. };
        if (lla_to_local(self, latitude, longitude, altitude, local) != 0)
                return NULL;
        return Py_BuildValue("(d,d,d)", local[0], local[1], local[2]);
}

/* Convert UTM coordinates to a cartesian position in local frame. */
static PyObject * topography_utm_to_local(
    TopographyObject * self, PyObject * args)
{
        double x, y, altitude;
        if (!PyArg_ParseTuple(args, "ddd", &x, &y, &altitude)) return NULL;
        double latitude = 0., longitude = 0.;
        if (turtle_projection_unproject(self->projection, x, y, &latitude,
                &longitude) != TURTLE_RETURN_SUCCESS)
                return NULL;
        double local[3] = { 0., 0., 0. };
        if (lla_to_local(self, latitude, longitude, altitude, local) != 0)
                return NULL;
        return Py_BuildValue("(d,d,d)", local[0], local[1], local[2]);
}

/* Convert GRAND angular coordinates to a local direction. */
static PyObject * topography_angular_to_local(
    TopographyObject * self, PyObject * args)
{
        /* Parse the arguments. */
        PyObject * position_obj = NULL;
        double theta, phi;
        if (!PyArg_ParseTuple(args, "Odd", &position_obj, &theta, &phi))
                return NULL;
        double local[3];
        if (parse_vector(position_obj, 3, local) != 0) return NULL;

        /* Compute the geodetic coordinates. */
        double latitude = 0., longitude = 0., altitude = 0.;
        if (local_to_lla(self, local, &latitude, &longitude, &altitude) != 0)
                return NULL;

        /* Compute the local direction. */
        double ecef[3];
        const double azimuth = -phi;
        const double elevation = 90. - theta;
        if (turtle_datum_direction(self->datum, latitude, longitude, azimuth,
                elevation, ecef) != TURTLE_RETURN_SUCCESS)
                return NULL;
        double direction[3];
        ecef_to_local(self, ecef, direction, 1);

        return Py_BuildValue(
            "(d,d,d)", direction[0], direction[1], direction[2]);
}

/* Encapsulate TURTLE calls to elevation in order to check for a flat
 * topography.
 */
static enum turtle_return ground_elevation(
    TopographyObject * self, double latitude, double longitude, double * z)
{
        if (self->flat) {
                *z = 0.;
                if ((fabs(latitude - self->latitude) > self->flat_size) ||
                    (fabs(longitude - self->longitude) > self->flat_size)) {
                        if (!self->catch)
                                handle_turtle_error(TURTLE_RETURN_PATH_ERROR,
                                    (turtle_caller_t *)turtle_datum_elevation);
                        return TURTLE_RETURN_PATH_ERROR;
                } else {
                        return TURTLE_RETURN_SUCCESS;
                }
        } else {
                return turtle_datum_elevation(
                    self->datum, latitude, longitude, z);
        }
}

/* Estimate the ground elevation in the local frame. */
static enum turtle_return ground_elevation_local(
    TopographyObject * self, double x, double y, double * zg)
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
                double latitude, longitude, altitude;
                if ((rc = local_to_lla(self, local, &latitude, &longitude,
                         &altitude)) != TURTLE_RETURN_SUCCESS)
                        return rc;
                double z;
                if ((rc = ground_elevation(self, latitude, longitude, &z)) !=
                    TURTLE_RETURN_SUCCESS)
                        return rc;

                /* Compute the corresponding local coordinates. */
                if ((rc = lla_to_local(self, latitude, longitude, z, new)) !=
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

/* Get the ground altitude in local frame coordinates or in geodetic ones. */
static PyObject * topography_ground_altitude(
    TopographyObject * self, PyObject * args, PyObject * kwargs)
{
        /* Parse the arguments. */
        double x, y;
        int geodetic = 0;

        static char * kwlist[] = { "x", "y", "geodetic", NULL };
        if (!PyArg_ParseTupleAndKeywords(
                args, kwargs, "dd|b", kwlist, &x, &y, &geodetic))
                return NULL;

        double z;
        if (geodetic) {
                if (ground_elevation(self, x, y, &z) != TURTLE_RETURN_SUCCESS)
                        return NULL;
        } else {
                if (ground_elevation_local(self, x, y, &z) !=
                    TURTLE_RETURN_SUCCESS)
                        return NULL;
        }

        return Py_BuildValue("d", z);
}

/* Compute the vertical distance to the ground. */
static enum turtle_return vertical_distance(
    TopographyObject * self, double * local, double * distance)
{
        enum turtle_return rc;
        double latitude, longitude, altitude;
        if ((rc = local_to_lla(self, local, &latitude, &longitude,
                 &altitude)) != TURTLE_RETURN_SUCCESS)
                return rc;
        double z;
        if ((rc = ground_elevation(self, latitude, longitude, &z)) !=
            TURTLE_RETURN_SUCCESS)
                return rc;

        *distance = altitude - z;
        return TURTLE_RETURN_SUCCESS;
}

/* Check if the given local position is above the ground. */
static enum turtle_return is_above(
    TopographyObject * self, double * local, int * status)
{
        enum turtle_return rc;
        double distance;
        if ((rc = vertical_distance(self, local, &distance)) !=
            TURTLE_RETURN_SUCCESS)
                return rc;
        *status = (distance > 0.);
        return TURTLE_RETURN_SUCCESS;
}

static PyObject * topography_is_above(
    TopographyObject * self, PyObject * position)
{
        double local[3];
        if (parse_vector(position, 3, local) != 0) return NULL;
        int status;
        if (is_above(self, local, &status) != 0) return NULL;
        return Py_BuildValue("b", status);
}

/* Compute the distance to the topography for the given point, or segment. */
static PyObject * topography_distance(
    TopographyObject * self, PyObject * args, PyObject * kwargs)
{
        /* Parse the arguments. */
        PyObject *position_obj, *direction_obj = NULL;
        double limit = -1.;

        static char * kwlist[] = { "position", "direction", "limit", NULL };
        if (!PyArg_ParseTupleAndKeywords(args, kwargs, "O|Od", kwlist,
                &position_obj, &direction_obj, &limit))
                return NULL;

        double position[3];
        if (parse_vector(position_obj, 3, position) != 0) return NULL;

        if (direction_obj == NULL) {
                /* Let us return the vertical distance to the ground. */
                double distance;
                if (vertical_distance(self, position, &distance) !=
                    TURTLE_RETURN_SUCCESS)
                        return NULL;
                return Py_BuildValue("d", distance);
        }

        /* Let us step along the given direction until the ground is
         * reached.
         */
        double direction[3];
        if (parse_vector(direction_obj, 3, direction) != 0) return NULL;

        const double resolution = 1E-02;
        int above;
        if (is_above(self, position, &above) != TURTLE_RETURN_SUCCESS)
                return NULL;
        double r0[3] = { position[0], position[1], position[2] };
        double distance = 0.;
        turtle_handler(NULL);
        self->catch = 1;
        for (;;) {
                double step;
                if (vertical_distance(self, r0, &step) != TURTLE_RETURN_SUCCESS)
                        break;
                step = 0.5 * fabs(step);
                if (step < resolution) step = resolution;
                double r1[3] = { r0[0] + step * direction[0],
                        r0[1] + step * direction[1],
                        r0[2] + step * direction[2] };
                int above1;
                if (is_above(self, r1, &above1) != TURTLE_RETURN_SUCCESS) break;
                if (above1 != above) {
                        /* A medium change was found. Let us refine it with a
                         * binary search.
                         */
                        while (step > resolution) {
                                double r2[3] = { 0.5 * (r0[0] + r1[0]),
                                        0.5 * (r0[1] + r1[1]),
                                        0.5 * (r0[2] + r1[2]) };
                                int above2 = 0;
                                is_above(self, r2, &above2);
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
                        turtle_handler(&handle_turtle_error);
                        self->catch = 0;
                        return Py_BuildValue("d", distance);
                } else if (fabs(r1[2] > 1E+04))
                        break;
                if (limit > 0.) {
                        distance += step;
                        if (distance > limit) break;
                }
                memcpy(r0, r1, sizeof(r0));
        }

        turtle_handler(&handle_turtle_error);
        self->catch = 0;
        Py_RETURN_NONE;
}

/* Compute the ground normal in local frame coordinates. */
static enum turtle_return ground_normal_local(
    TopographyObject * self, double x, double y, double * n)
{
        /* Estimate the local slope using an adaptation of S. Le Coz's
         * algorithm.
         */
        enum turtle_return rc;
        const double ds = 10.;
        double zx0, zx1, zy0, zy1;
        if ((rc = ground_elevation_local(self, x - 0.5 * ds, y, &zx0)) !=
            TURTLE_RETURN_SUCCESS)
                return rc;
        if ((rc = ground_elevation_local(self, x + 0.5 * ds, y, &zx1)) !=
            TURTLE_RETURN_SUCCESS)
                return rc;
        if ((rc = ground_elevation_local(self, x, y - 0.5 * ds, &zy0)) !=
            TURTLE_RETURN_SUCCESS)
                return rc;
        if ((rc = ground_elevation_local(self, x, y + 0.5 * ds, &zy1)) !=
            TURTLE_RETURN_SUCCESS)
                return rc;

        double u[3] = { ds, 0., zx1 - zx0 };
        double v[3] = { 0., ds, zy1 - zy0 };
        n[0] = u[1] * v[2] - u[2] * v[1];
        n[1] = u[2] * v[0] - u[0] * v[2];
        n[2] = u[0] * v[1] - u[1] * v[0];
        double d = n[0] * n[0] + n[1] * n[1] + n[2] * n[2];
        if (d <= FLT_EPSILON)
                d = 0.;
        else
                d = 1. / sqrt(d);
        n[0] *= d;
        n[1] *= d;
        n[2] *= d;

        return TURTLE_RETURN_SUCCESS;
}

/* Compute the normal to the ground, and optionnaly the corresponding angles. */
static PyObject * topography_ground_normal(
    TopographyObject * self, PyObject * args, PyObject * kwargs)
{
        /* Parse the arguments. */
        double x, y;
        int geodetic = 0, angles = 0;

        static char * kwlist[] = { "x", "y", "geodetic", "angles", NULL };
        if (!PyArg_ParseTupleAndKeywords(
                args, kwargs, "dd|bb", kwlist, &x, &y, &geodetic, &angles))
                return NULL;

        double n[3] = { 0., 0., 0. };
        if (self->flat) {
                /* Compute the local verticale. */
                if (!geodetic) {
                        double z;
                        if (ground_elevation_local(self, x, y, &z) !=
                            TURTLE_RETURN_SUCCESS)
                                return NULL;

                        double local[3] = { x, y, z }, tmp;
                        local_to_lla(self, local, &x, &y, &tmp);
                }

                double ecef[3];
                if (turtle_datum_direction(self->datum, x, y, 0., 90., ecef) !=
                    TURTLE_RETURN_SUCCESS)
                        return NULL;
                ecef_to_local(self, ecef, n, 1);
        } else {
                /* Estimate the local slope of the ground. */
                if (geodetic) {
                        double z;
                        if (ground_elevation(self, x, y, &z) !=
                            TURTLE_RETURN_SUCCESS)
                                return NULL;

                        double local[3];
                        if (lla_to_local(self, x, y, z, local) !=
                            TURTLE_RETURN_SUCCESS)
                                return NULL;
                        x = local[0];
                        y = local[1];
                }

                if (ground_normal_local(self, x, y, n) != TURTLE_RETURN_SUCCESS)
                        return NULL;
        }

        if (angles) {
                const double phi = atan2(n[1], n[0]) * 180. / M_PI;
                const double theta = acos(n[2]) * 180. / M_PI;
                return Py_BuildValue(
                    "(d, d, d), d, d", n[0], n[1], n[2], theta, phi);
        } else {
                return Py_BuildValue("(d, d, d)", n[0], n[1], n[2]);
        }
}

/* Register the Topography methods. */
static PyMethodDef topography_methods[] = {
        { "local_to_ecef", (PyCFunction)topography_local_to_ecef, METH_VARARGS,
            "Convert a cartesian position in local frame to an ECEF one" },
        { "local_to_lla", (PyCFunction)topography_local_to_lla, METH_O,
            "Convert a cartesian position in local frame to a geodetic one" },
        { "local_to_utm", (PyCFunction)topography_local_to_utm, METH_O,
            "Convert a cartesian position in local frame to UTM coordinates" },
        { "local_to_angular", (PyCFunction)topography_local_to_angular,
            METH_VARARGS, "Convert a local direction to angular coordinates" },
        { "ecef_to_local", (PyCFunction)topography_ecef_to_local, METH_VARARGS,
            "Convert a cartesian position in ECEF frame to a local one" },
        { "lla_to_local", (PyCFunction)topography_lla_to_local, METH_VARARGS,
            "Convert a geodetic position to a cartesian one in local frame" },
        { "utm_to_local", (PyCFunction)topography_utm_to_local, METH_VARARGS,
            "Convert UTM coordinates to a cartesian position in local frame" },
        { "angular_to_local", (PyCFunction)topography_angular_to_local,
            METH_VARARGS, "Convert angular coordinates to a local direction" },
        { "ground_altitude", (PyCFunction)topography_ground_altitude,
            METH_KEYWORDS,
            "Get the ground altitude in local frame coordinates "
            "or in geodetic ones" },
        { "ground_normal", (PyCFunction)topography_ground_normal, METH_KEYWORDS,
            "Get the normal to the ground in local frame coordinates." },
        { "is_above", (PyCFunction)topography_is_above, METH_O,
            "Check if the given local position is above the ground" },
        { "distance", (PyCFunction)topography_distance, METH_KEYWORDS,
            "Compute the distance to the topography for the given point, "
            "or segment" },
        { NULL }
};

/* The topography type. */
static PyTypeObject TopographyType = { PyVarObject_HEAD_INIT(
                                           NULL, 0) "grand_tour.Topography",
        sizeof(TopographyObject), .tp_flags = Py_TPFLAGS_DEFAULT,
        .tp_doc = "Topography provider in a cartesian local frame",
        .tp_init = (initproc)topography_init,
        .tp_dealloc = (destructor)topography_destroy,
        .tp_methods = topography_methods };

/* Initialise the module. */
PyMODINIT_FUNC initgrand_tour(void)
{
        /* Initialise the Topography type. */
        TopographyType.tp_new = PyType_GenericNew;
        if (PyType_Ready(&TopographyType) < 0) return;

        /* Initialise the module. */
        PyObject * module;
        if ((module = Py_InitModule("grand_tour", NULL)) == NULL) return;

        /* Register the topography type. */
        Py_INCREF(&TopographyType);
        PyModule_AddObject(module, "Topography", (PyObject *)&TopographyType);

        /* Register the TurtleError exception. */
        TurtleError = PyErr_NewException("grand_tour.TurtleError", NULL, NULL);
        Py_INCREF(TurtleError);
        PyModule_AddObject(module, "TurtleError", TurtleError);

        /* Initialise the turtle library. */
        turtle_initialise(&handle_turtle_error);
}
