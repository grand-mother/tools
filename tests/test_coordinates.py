# -*- coding: utf-8 -*-
"""
Unit tests for the tools.coordinates module
"""

import unittest

from tools.coordinates import ENU, GeodeticRepresentation

import numpy
import astropy.units as u
from astropy.coordinates import CartesianRepresentation, EarthLocation, ITRS,  \
                                SkyCoord


class CoordinatesTest(unittest.TestCase):
    """Unit tests for the coordinates module"""


    def assertQuantity(self, x, y, unit, tol):
        """Check that two astropy.Quantities are consistent"""
        self.assertAlmostEqual((x / unit).value, (y / unit).value, tol)


    def test_geodetic(self):
        settings = {"lat" : 45.0 * u.deg, "lon" : 3.0 * u.deg,
                    "height" : 1000. * u.m}
        location = EarthLocation(**settings)
        cart = location.itrs.cartesian

        # Check the Cartesian generator
        r = GeodeticRepresentation.from_cartesian(cart)

        self.assertQuantity(r.latitude, settings["lat"], u.deg, 11)
        self.assertQuantity(r.longitude, settings["lon"], u.deg, 11)
        self.assertQuantity(r.height, settings["height"], u.m, 6)

        # Check the Cartesian conversion
        r = GeodeticRepresentation(settings["lat"], settings["lon"],
                                   settings["height"]).to_cartesian()

        self.assertQuantity(r.x, cart.x, u.m, 6)
        self.assertQuantity(r.y, cart.y, u.m, 6)
        self.assertQuantity(r.z, cart.z, u.m, 6)

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
            self.assertQuantity(r.latitude[i], settings["lat"], u.deg, 11)
            self.assertQuantity(r.longitude[i], settings["lon"], u.deg, 11)
            self.assertQuantity(r.height[i], settings["height"], u.m, 6)

        # Check the vetcorized Cartesian conversion
        r = GeodeticRepresentation(latitude, longitude, height).to_cartesian()

        for i in range(n):
            self.assertQuantity(r.x[i], cart.x[i], u.m, 6)
            self.assertQuantity(r.y[i], cart.y[i], u.m, 6)
            self.assertQuantity(r.z[i], cart.z[i], u.m, 6)


    def test_enu(self):
        settings = {"lat" : 45.0 * u.deg, "lon" : 3.0 * u.deg,
                    "height" : 1000. * u.m}
        location = EarthLocation(**settings)
        itrs = location.itrs

        # Check the constructor & to ITRS transform
        enu = ENU(x=0 * u.m, y=0 * u.m, z=0 * u.m, location=location) 
        r = enu.transform_to(ITRS)

        self.assertQuantity(r.x, itrs.x, u.m, 6)
        self.assertQuantity(r.y, itrs.y, u.m, 6)
        self.assertQuantity(r.z, itrs.z, u.m, 6)

        # Check the from ITRS transform
        enu = ENU(location=location)
        enu = itrs.transform_to(enu)

        self.assertQuantity(enu.x, 0, u.m, 6)
        self.assertQuantity(enu.y, 0, u.m, 6)
        self.assertQuantity(enu.z, 0, u.m, 6)

        # Check the affine transform
        points = ((0, 0, 1), (1, 0, 0), (0, 1, 0), (1, 1, 0), (1, 1, 1),
                  (0, 1, 1))
        for point in (points):
            cart = CartesianRepresentation(x=point[1], y=point[0], z=point[2],
                                           unit=u.m)
            altaz = SkyCoord(cart, frame="altaz", location=location,
                             obstime="2019-01-01")
            itrs0 = altaz.transform_to(ITRS)

            cart = CartesianRepresentation(x=point[0], y=point[1], z=point[2],
                                           unit=u.m)
            enu = ENU(cart, location=location)
            itrs1 = enu.transform_to(ITRS)

            self.assertQuantity(itrs0.x, itrs1.x, u.m, 4)
            self.assertQuantity(itrs0.y, itrs1.y, u.m, 4)
            self.assertQuantity(itrs0.z, itrs1.z, u.m, 4)

        # Check the orientation
        point = (1, -1, 2)
        cart = CartesianRepresentation(x=point[0], y=point[1], z=point[2],
                                       unit=u.m)
        altaz = SkyCoord(cart, frame="altaz", location=location,
                         obstime="2019-01-01")
        itrs0 = altaz.transform_to(ITRS)

        for (orientation, sign) in ((("N", "E", "U"), (1, 1, 1)),
                                    (("N", "E", "D"), (1, 1, -1)),
                                    (("S", "E", "U"), (-1, 1, 1)),
                                    (("N", "W", "U"), (1, -1, 1))):
            cart = CartesianRepresentation(x=sign[0] * point[0],
                y=sign[1] * point[1], z=sign[2] * point[2], unit=u.m)
            enu = ENU(cart, location=location, orientation=orientation)
            itrs1 = enu.transform_to(ITRS)

            self.assertQuantity(itrs0.x, itrs1.x, u.m, 4)
            self.assertQuantity(itrs0.y, itrs1.y, u.m, 4)
            self.assertQuantity(itrs0.z, itrs1.z, u.m, 4)


if __name__ == "__main__":
    unittest.main()
