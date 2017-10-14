# -*- coding: utf-8 -*-
#
#  Copyright (C) 2017 Université Clermont Auvergne, CNRS/IN2P3, LPC
#  Author: Valentin NIESS (niess@in2p3.fr)
#
#  GRAND TOpography with tURtle (GRAND-TOUR)
#
#  This program is free software: you can redistribute it and/or modify
#  it under the terms of the GNU Lesser General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU Lesser General Public License for more details.
#
#  You should have received a copy of the GNU Lesser General Public License
#  along with this program.  If not, see <http://www.gnu.org/licenses/>

from ctypes import *

def _initialise_turtle():
    # Load the library and initialise it.
    turtle = cdll.LoadLibrary("libturtle.so")
    turtle.turtle_initialise(0)

    # Prototype of the `turtle_strerror` library function.
    turtle.turtle_strerror.argtypes = (c_int,)
    turtle.turtle_strerror.restype = c_char_p

    # Prototype of the `turtle_datum_create` library function.
    turtle.turtle_datum_create.argtypes = (c_char_p, c_int, c_void_p,
      c_void_p, POINTER(c_void_p))

    # Prototype of the `turtle_datum_destroy` library function.
    turtle.turtle_datum_destroy.argtypes = (POINTER(c_void_p),)

    # Prototype of the `turtle_datum_elevation` library function.
    turtle.turtle_datum_elevation.argtypes = (c_void_p, c_double, c_double,
      POINTER(c_double))

    # Prototype of the `turtle_datum_ecef` library function.
    turtle.turtle_datum_ecef.argtypes = (c_void_p, c_double, c_double,
      c_double, POINTER(c_double))

    # Prototype of the `turtle_datum_geodetic` library function.
    turtle.turtle_datum_geodetic.argtypes = (c_void_p, POINTER(c_double),
      POINTER(c_double), POINTER(c_double), POINTER(c_double))

    # Prototype of the `turtle_datum_direction` library function.
    turtle.turtle_datum_direction.argtypes = (c_void_p, c_double, c_double,
      c_double, c_double, POINTER(c_double))

    # Prototype of the `turtle_datum_horizontal` library function.
    turtle.turtle_datum_horizontal.argtypes = (c_void_p, c_double, c_double,
       POINTER(c_double),  POINTER(c_double), POINTER(c_double))

    # Prototype of the `turtle_projection_create` library function.
    turtle.turtle_projection_create.argtypes = (c_char_p, POINTER(c_void_p))

    # Prototype of the `turtle_projection_destroy` library function.
    turtle.turtle_projection_destroy.argtypes = (POINTER(c_void_p),)

    # Prototype of the `turtle_projection_project` library function.
    turtle.turtle_projection_project.argtypes = (c_void_p, c_double, c_double,
      POINTER(c_double), POINTER(c_double))

    # Prototype of the `turtle_projection_unproject` library function.
    turtle.turtle_projection_unproject.argtypes = (c_void_p, c_double, c_double,
      POINTER(c_double), POINTER(c_double))

    return turtle

_turtle = _initialise_turtle()

class TurtleError(Exception):
    """Custom exception for a TURTLE library error.
    """
    def __init__(self, rc):
        message = _turtle.turtle_strerror(rc)
        super(TurtleError, self).__init__(message)
        self.rc = rc

class Datum(object):
    """Encapsulation of a TURTLE datum.
    """
    def __init__(self, path, stack_size=4):
        null = c_void_p(0)
        self._datum = c_void_p(0)
        rc = _turtle.turtle_datum_create(
          path, stack_size, null, null, byref(self._datum))
        if rc != 0: raise TurtleError(rc)

    def __del__(self):
        if self._datum is None: return
        _turtle.turtle_datum_destroy(byref(self._datum))
        self._datum = None

    def ground_altitude(self, latitude, longitude):
        """Get the ground altitude at a given geodetic location.
        """
        z = c_double(0.)
        rc = _turtle.turtle_datum_elevation(
          self._datum, latitude, longitude, byref(z))
        if rc != 0: raise TurtleError(rc)
        return z.value

    def lla_to_ecef(self, latitude, longitude, altitude):
        """Convert a geodetic position to a cartesian one, in ECEF.
        """
        r = (3 * c_double)(0., 0., 0.)
        rc = _turtle.turtle_datum_ecef(
          self._datum, latitude, longitude, altitude, r)
        if rc != 0: raise TurtleError(rc)
        return tuple(x for x in r)

    def ecef_to_lla(self, position):
        """Convert a cartesian position in ECEF to a geodetic one.
        """
        r = (3 * c_double)(*position)
        latitude, longitude, altitude = c_double(0), c_double(0), c_double(0)
        rc = _turtle.turtle_datum_geodetic(
          self._datum, r, byref(latitude), byref(longitude), byref(altitude))
        if rc != 0: raise TurtleError(rc)
        return (latitude.value, longitude.value, altitude.value)

    def horizontal_to_ecef(self, latitude, longitude, azimuth, elevation):
        """Convert an horizontal direction to an ECEF one.
        """
        r = (3 * c_double)(0., 0., 0.)
        rc = _turtle.turtle_datum_direction(
          self._datum, latitude, longitude, azimuth, elevation, r)
        if rc != 0: raise TurtleError(rc)
        return tuple(x for x in r)

    def ecef_to_horizontal(self, latitude, longitude, position):
        """Convert an ECEF direction to an horizontal one.
        """
        r = (3 * c_double)(*position)
        azimuth, elevation = c_double(0.), c_double(0.)
        rc = _turtle.turtle_datum_horizontal(
          self._datum, latitude, longitude, r, byref(azimuth), byref(elevation))
        if rc != 0: raise TurtleError(rc)
        return azimuth.value, elevation.value

class Projection(object):
    """Encapsulation of a TURTLE projection.
    """
    def __init__(self, name):
        self._projection = c_void_p(0)
        rc = _turtle.turtle_projection_create(name, byref(self._projection))
        if rc != 0: raise TurtleError(rc)

    def __del__(self):
        if self._projection is None: return
        _turtle.turtle_projection_destroy(byref(self._projection))
        self._projection = None

    def project(self, latitude, longitude):
        """Convert a geodetic position to map coordinates.
        """
        x, y = c_double(0.), c_double(0.)
        rc = _turtle.turtle_projection_project(self._projection,
          latitude, longitude, byref(x), byref(y))
        if rc != 0: raise TurtleError(rc)
        return x.value, y.value

    def unproject(self, x, y):
        """Convert map coordinates to a geodetic position.
        """
        latitude, longitude = c_double(0.), c_double(0.)
        rc = _turtle.turtle_projection_unproject(self._projection,
          x, y, byref(latitude), byref(longitude))
        if rc != 0: raise TurtleError(rc)
        return latitude.value, longitude.value

class Topography(Datum, Projection):
    """Topography provider in a cartesian local frame.
    """
    def __init__(self, latitude, longitude, path, stack_size=4):
        Datum.__init__(self, path, stack_size)

        # Set the local frame.
        self.origin = self.lla_to_ecef(latitude, longitude, 0.)
        self.base = [self.horizontal_to_ecef(latitude, longitude, 90., 0.),
          self.horizontal_to_ecef(latitude, longitude, 0., 0.),
          self.horizontal_to_ecef(latitude, longitude, 0., 90.)]

        # Set the UTM projection.
        zone = (int((longitude + 180) / 6) % 60) + 1
        if latitude >= 0.: hemisphere = "N"
        else: hemisphere = "S"
        Projection.__init__(self, "UTM {:}{:}".format(zone, hemisphere))

    def __del__(self):
        Projection.__del__(self)
        Datum.__del__(self)

    def local_to_lla(self, position):
        """Convert a cartesian position in local frame to a geodetic one.
        """
        r = [sum(self.base[j][i] * position[j] for j in xrange(3)) +
             self.origin[i] for i in xrange(3)]
        return self.ecef_to_lla(r)

    def local_to_utm(self, position):
        """Convert a cartesian position in local frame to UTM coordinates.
        """
        lla = self.local_to_lla(position)
        x, y = self.project(lla[0], lla[1])
        return x, y, lla[2]

    def lla_to_local(self, latitude, longitude, altitude):
        """Convert a geodetic position to a cartesian one in local frame.
        """
        r = self.lla_to_ecef(latitude, longitude, altitude)
        return [sum(self.base[i][j] * (r[j] - self.origin[j])
                for j in xrange(3)) for i in xrange(3)]

    def utm_to_local(self, x, y, altitude):
        """Convert UTM coordinates to a cartesian position in local frame.
        """
        latitude, longitude = self.unproject(x, y)
        return self.lla_to_local(latitude, longitude, altitude)

    def ground_altitude(self, x, y, geodetic=False):
        """Get the ground altitude in local frame coordinates or geodetic ones.
        """
        if geodetic:
            return super(Topography, self).ground_altitude(x, y)

        # Compute the altitude in local coordinates using an iterative method.
        # We need to find the ground position that projects at (x, y). Let us
        # Initialise the search on the ellipsoid, i.e. zero altitude.
        local = [x, y, 0.]

        for _ in xrange(5):
            # Get the ground altitude for the current guess.
            lla = self.local_to_lla(local)
            ground = super(Topography, self).ground_altitude(lla[0], lla[1])

            # Compute the corresponding local ccordinates.
            new = self.lla_to_local(lla[0], lla[1], ground)

            # Check for convergence and update.
            dx, dy = x - new[0], y - new[1]
            if (abs(dx) < 1E+03) and (abs(dy) < 1E-03): break
            local[0] += dx
            local[1] += dy
        return new[2]

    def is_above(self, position):
        """Check if the given local position is above the ground.
        """
        return self._vertical_distance(position) > 0.

    def distance(self, position, direction=None, limit=None):
        """Compute the distance to the topography for the given point,
           or segment.
        """
        if direction:
            # Let us step along the given direction until the ground is
            # reached.
            above = self.is_above(position)
            r0 = [c for c in position]
            distance = 0.
            while True:
                step = 0.5 * abs(self._vertical_distance(position))
                if step < 1E-02: step = 1E-02
                r1 = [r0[i] + step * direction[i] for i in xrange(3)]
                try: above1 = self.is_above(r1)
                except TurtleError: return None
                if above1 != above:
                    # A medium change was found. Let us refine with a
                    # binary search.
                    while step > 1E-02:
                        r2 = [0.5 * (r0[i] + r1[i]) for i in xrange(3)]
                        if self.is_above(r2) == above: r0 = r2
                        else: r1 = r2
                        step *= 0.5
                    distance = sum((r1[i] - position[i])**2
                      for i in xrange(3))**0.5
                    if limit and (distance > limit): return None
                    return distance
                elif abs(r1[2]) > 1E+04: return None
                if limit:
                    distance += step
                    if distance > limit: return None
                r0 = r1
        else:
            # Let us return the vertical distance to the ground.
            return self._vertical_distance(position)

    def _vertical_distance(self, position):
        """Compute the vertical distance to the ground."""
        lla = self.local_to_lla(position)
        ground = super(Topography, self).ground_altitude(lla[0], lla[1])
        return lla[2] - ground
