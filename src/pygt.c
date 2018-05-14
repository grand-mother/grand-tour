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
#include "grand-tour.h"

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

            /* Low level toopography object */
            struct gt_topography topography;
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
        self->topography.catch = 0;

        const int rc = gt_initialise(&self->topography, latitude, longitude,
            path, stack_size, &handle_turtle_error, NULL);
        if (rc < 0) {
                PyErr_SetString(
                    PyExc_ValueError, "invalid size for flat topography");
                return -1;
        } else if (rc != TURTLE_RETURN_SUCCESS)
                return -1;
        return 0;
}

/* Free a topography object. */
static void topography_destroy(TopographyObject * self)
{
        turtle_datum_destroy(&self->topography.datum);
        turtle_projection_destroy(&self->topography.projection);
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
        gt_to_ecef(&self->topography, local, vector, ecef);

        return Py_BuildValue("(d,d,d)", ecef[0], ecef[1], ecef[2]);
}

static PyObject * topography_local_to_lla(
    TopographyObject * self, PyObject * position)
{
        double local[3];
        if (parse_vector(position, 3, local) != 0) return NULL;
        double lla[3] = { 0., 0., 0. };
        if (gt_to_lla(&self->topography, local, lla) != TURTLE_RETURN_SUCCESS)
                return NULL;
        return Py_BuildValue("ddd", lla[0], lla[1], lla[2]);
}

/* Convert a cartesian position in local frame to UTM coordinates. */
static PyObject * topography_local_to_utm(
    TopographyObject * self, PyObject * position)
{
        double local[3];
        if (parse_vector(position, 3, local) != 0) return NULL;
        double utm[3];
        if (gt_to_utm(&self->topography, local, utm) != TURTLE_RETURN_SUCCESS)
                return NULL;
        return Py_BuildValue("ddd", utm[0], utm[1], utm[2]);
}

/* Convert a local direction to GRAND angular coordinates. */
static PyObject * topography_local_to_angular(
    TopographyObject * self, PyObject * args)
{
        /* Parse the arguments. */
        PyObject *position_obj = NULL, *direction_obj = NULL;
        if (!PyArg_ParseTuple(args, "OO", &position_obj, &direction_obj))
                return NULL;
        double position[3];
        if (parse_vector(position_obj, 3, position) != 0) return NULL;
        double direction[3];
        if (parse_vector(direction_obj, 3, direction) != 0) return NULL;

        double angular[2] = { 0., 0. };
        if (gt_to_angular(&self->topography, position, direction, angular) !=
            TURTLE_RETURN_SUCCESS)
                return NULL;

        return Py_BuildValue("dd", angular[0], angular[1]);
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
        gt_from_ecef(&self->topography, ecef, vector, local);

        return Py_BuildValue("(d,d,d)", local[0], local[1], local[2]);
}

static PyObject * topography_lla_to_local(
    TopographyObject * self, PyObject * args)
{
        double lla[3];
        if (!PyArg_ParseTuple(args, "ddd", lla, lla + 1, lla + 2)) return NULL;
        double local[3] = { 0., 0., 0. };
        if (gt_from_lla(&self->topography, lla, local) != 0) return NULL;
        return Py_BuildValue("(d,d,d)", local[0], local[1], local[2]);
}

/* Convert UTM coordinates to a cartesian position in local frame. */
static PyObject * topography_utm_to_local(
    TopographyObject * self, PyObject * args)
{
        double utm[3] = { 0., 0., 0. };
        if (!PyArg_ParseTuple(args, "ddd", utm, utm + 1, utm + 2)) return NULL;

        int rc;
        double local[3] = { 0., 0., 0. };
        if ((rc = gt_from_utm(&self->topography, utm, local)) !=
            TURTLE_RETURN_SUCCESS)
                return NULL;

        return Py_BuildValue("(d,d,d)", local[0], local[1], local[2]);
}

/* Convert GRAND angular coordinates to a local direction. */
static PyObject * topography_angular_to_local(
    TopographyObject * self, PyObject * args)
{
        /* Parse the arguments. */
        PyObject * position_obj = NULL;
        double angular[2];
        if (!PyArg_ParseTuple(args, "Odd", &position_obj, angular, angular + 1))
                return NULL;
        double position[3];
        if (parse_vector(position_obj, 3, position) != 0) return NULL;

        double direction[3] = { 0., 0., 0. };
        if (gt_from_angular(&self->topography, position, angular, direction) !=
            TURTLE_RETURN_SUCCESS)
                return NULL;

        return Py_BuildValue(
            "(d,d,d)", direction[0], direction[1], direction[2]);
}

/* Get the ground altitude in local frame coordinates or in geodetic ones. */
static PyObject * topography_ground_altitude(
    TopographyObject * self, PyObject * args, PyObject * kwargs)
{
        /* Parse the arguments. */
        double position[2];
        int geodetic = 0;

        static char * kwlist[] = { "x", "y", "geodetic", NULL };
        if (!PyArg_ParseTupleAndKeywords(args, kwargs, "dd|b", kwlist, position,
                position + 1, &geodetic))
                return NULL;

        double altitude;
        int rc;
        if ((rc = gt_ground_altitude(&self->topography, position, geodetic,
                 &altitude)) != TURTLE_RETURN_SUCCESS)
                return NULL;

        return Py_BuildValue("d", altitude);
}

static PyObject * topography_is_above(
    TopographyObject * self, PyObject * position)
{
        double local[3];
        if (parse_vector(position, 3, local) != 0) return NULL;
        int status;
        if (gt_ground_above(&self->topography, local, &status) !=
            TURTLE_RETURN_SUCCESS)
                return NULL;
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

        double *direction = NULL, tmp[3];
        if (direction_obj != NULL) {
                if (parse_vector(direction_obj, 3, tmp) != 0) return NULL;
                direction = tmp;
        }

        double distance;
        if (gt_ground_distance(&self->topography, position, direction, limit,
                &distance) != TURTLE_RETURN_SUCCESS)
                return NULL;

        if (distance >= 0.) return Py_BuildValue("d", distance);
        Py_RETURN_NONE;
}

/* Compute the normal to the ground, and optionnaly the corresponding angles. */
static PyObject * topography_ground_normal(
    TopographyObject * self, PyObject * args, PyObject * kwargs)
{
        /* Parse the arguments. */
        double position[2], step = -1.;
        int geodetic = 0, compute_angles = 0;

        static char * kwlist[] = { "x", "y", "geodetic", "step", "angles",
                NULL };
        if (!PyArg_ParseTupleAndKeywords(args, kwargs, "dd|bdb", kwlist,
                position, position + 1, &geodetic, &step, &compute_angles))
                return NULL;

        double normal[3] = { 0., 0., 0. }, angles_[2] = { 0., 0. };
        double * angles_v = compute_angles ? angles_ : NULL;
        if (gt_ground_normal(&self->topography, position, geodetic, step,
                normal, angles_v) != TURTLE_RETURN_SUCCESS)
                return NULL;

        if (compute_angles) {
                return Py_BuildValue("(d, d, d), d, d", normal[0], normal[1],
                    normal[2], angles_v[0], angles_v[1]);
        } else {
                return Py_BuildValue(
                    "(d, d, d)", normal[0], normal[1], normal[2]);
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
