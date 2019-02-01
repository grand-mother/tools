# -*- coding: utf-8 -*-
"""
Unit tests for the tools.coordinates module
"""

import unittest

from tools.coordinates import GeodeticRepresentation

import numpy
import astropy.units as u
from astropy.coordinates import CartesianRepresentation, EarthLocation, ITRS


class CoordinatesTest(unittest.TestCase):
    """Unit tests for the coordinates module"""

    def test_geodetic(self):
        settings = {"lat" : 45.0 * u.deg, "lon" : 3.0 * u.deg,
                    "height" : 1000. * u.m}
        location = EarthLocation(**settings)
        cart = location.itrs.cartesian

        # Check the Cartesian generator
        r = GeodeticRepresentation.from_cartesian(cart)

        def assertQuantity(x, y, unit, tol):
            self.assertAlmostEqual((x / unit).value, (y / unit).value, tol)

        assertQuantity(r.latitude, settings["lat"], u.deg, 11)
        assertQuantity(r.longitude, settings["lon"], u.deg, 11)
        assertQuantity(r.height, settings["height"], u.m, 6)

        # Check the Cartesian conversion
        r = GeodeticRepresentation(settings["lat"], settings["lon"],
                                   settings["height"]).to_cartesian()

        assertQuantity(r.x, cart.x, u.m, 6)
        assertQuantity(r.y, cart.y, u.m, 6)
        assertQuantity(r.z, cart.z, u.m, 6)

        # Vectorize the test point
        n = 10
        vectorize = lambda v: u.Quantity(
            numpy.repeat((v[0] / v[1]).value, n), v[1])
        x, y, z = map(vectorize, ((cart.x, u.m), (cart.y, u.m), (cart.z, u.m)))
        cart = CartesianRepresentation(x=x, y=y, z=z)
        latitude, longitude, height = map(vectorize, ((settings["lat"], u.deg),
            (settings["lon"], u.deg), (settings["height"], u.m)))

        # Check the vectorized Cartesian generator
        r = GeodeticRepresentation.from_cartesian(cart)

        for i in range(n):
            assertQuantity(r.latitude[i], settings["lat"], u.deg, 11)
            assertQuantity(r.longitude[i], settings["lon"], u.deg, 11)
            assertQuantity(r.height[i], settings["height"], u.m, 6)

        # Check the vetcorized Cartesian conversion
        r = GeodeticRepresentation(latitude, longitude, height).to_cartesian()

        for i in range(n):
            assertQuantity(r.x[i], cart.x[i], u.m, 6)
            assertQuantity(r.y[i], cart.y[i], u.m, 6)
            assertQuantity(r.z[i], cart.z[i], u.m, 6)


if __name__ == "__main__":
    unittest.main()
