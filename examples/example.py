#!/usr/bin/env python
import time
from grand_tour import Topography
import numpy
import pylab

# Instanciate the topography.
topo = Topography(42.928056, 86.741667, "share/topography")

# Sample the local altitude over a mesh.
t0 = time.time()
print "o Sampling the local altitude ..."
x = y = numpy.linspace(-100E+03, 100E+03, 201)
z = numpy.zeros((len(y), len(x)))
for i, yi in enumerate(y):
    for j, xj in enumerate(x):
        z[i, j] = topo.ground_altitude(xj, yi)
print "  --> Done in {:.1f}s".format(time.time() - t0)

# Sample the local slope.
t0 = time.time()
print "o Sampling the slope vertical angle ..."
slope = numpy.zeros((len(y), len(x)))
for i, yi in enumerate(y):
    for j, xj in enumerate(x):
        slope[i, j] = topo.ground_normal(xj, yi, angles=True)[1]
print "  --> Done in {:.1f}s".format(time.time() - t0)

# Sample the geodetic altitude over a mesh in UTM coordinates.
t0 = time.time()
print "o Sampling the UTM altitude ..."
xUTM, yUTM, _ = topo.local_to_utm((0., 0., 0.))
zUTM = numpy.zeros((len(y), len(x)))
for i, yi in enumerate(y):
    for j, xj in enumerate(x):
        local = topo.utm_to_local(xj + xUTM, yi + yUTM, 0.)
        latitude, longitude, _ = topo.local_to_lla(local)
        zUTM[i, j] = topo.ground_altitude(latitude, longitude, geodetic=True)
print "  --> Done in {:.1f}s".format(time.time() - t0)

# Plot the result.
cmap = pylab.get_cmap("terrain")
opts = { "cmap" : cmap, "vmin" : 0., "vmax" : 5. }

pylab.figure()
pylab.pcolor(x * 1E-03, y * 1E-03, z * 1E-03, **opts)
pylab.colorbar()
pylab.xlabel("local x (km)")
pylab.ylabel("local y (km)")
pylab.title("local altitude (km)")

pylab.figure()
pylab.pcolor(x * 1E-03, y * 1E-03, slope, cmap=cmap, vmin=0., vmax=90.)
pylab.colorbar()
pylab.xlabel("local x (km)")
pylab.ylabel("local y (km)")
pylab.title("slope vertical angle (deg)")

pylab.figure()
pylab.pcolor(x * 1E-03, y * 1E-03, zUTM * 1E-03, **opts)
pylab.colorbar()
pylab.xlabel("UTM $\Delta$x (km)")
pylab.ylabel("UTM $\Delta$y (km)")
pylab.title("UTM altitude (km)")

pylab.show()
