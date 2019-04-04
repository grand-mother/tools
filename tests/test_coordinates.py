# -*- coding: utf-8 -*-
"""
Unit tests for the grand_tools.coordinates module
"""

import unittest

from grand_tools.coordinates import ECEF, ENU, GeodeticRepresentation,         \
                                    HorizontalRepresentation

import numpy
import astropy.units as u
from astropy.coordinates import CartesianRepresentation, EarthLocation, ITRS,  \
                                SkyCoord


class CoordinatesTest(unittest.TestCase):
    """Unit tests for the coordinates module"""

    def __init__(self, *args):
        super().__init__(*args)

        self.obstime = "2019-01-01"
        self.location = EarthLocation(lat=45.0 * u.deg, lon=3.0 * u.deg,
                                      height=1000. * u.m)

    def assertQuantity(self, x, y, unit, tol):
        """Check that two astropy.Quantities are consistent"""
        self.assertAlmostEqual((x / unit).value, (y / unit).value, tol)


    def assertCartesian(self, x, y, unit, tol):
        """Check that two astropy Cartesian coordinates are consistent"""
        self.assertQuantity(x.x, y.x, unit, tol)
        self.assertQuantity(x.y, y.y, unit, tol)
        self.assertQuantity(x.z, y.z, unit, tol)


    def test_geodetic(self):
        cart = self.location.itrs.cartesian

        # Check the Cartesian generator
        r = GeodeticRepresentation.from_cartesian(cart)

        self.assertQuantity(r.latitude, self.location.lat, u.deg, 11)
        self.assertQuantity(r.longitude, self.location.lon, u.deg, 11)
        self.assertQuantity(r.height, self.location.height, u.m, 6)

        # Check the Cartesian conversion
        r = GeodeticRepresentation(self.location.lat, self.location.lon,
                                   self.location.height).to_cartesian()

        self.assertQuantity(r.x, cart.x, u.m, 6)
        self.assertQuantity(r.y, cart.y, u.m, 6)
        self.assertQuantity(r.z, cart.z, u.m, 6)

        # Vectorize the test point
        n = 10
        vectorize = lambda v: u.Quantity(
            numpy.repeat((v[0] / v[1]).value, n), v[1])
        x, y, z = map(vectorize, ((cart.x, u.m), (cart.y, u.m), (cart.z, u.m)))
        cart = CartesianRepresentation(x=x, y=y, z=z)
        latitude, longitude, height = map(vectorize, (
            (self.location.lat, u.deg), (self.location.lon, u.deg),
            (self.location.height, u.m)))

        # Check the vectorized Cartesian generator
        r = GeodeticRepresentation.from_cartesian(cart)

        for i in range(n):
            self.assertQuantity(r.latitude[i], self.location.lat, u.deg, 11)
            self.assertQuantity(r.longitude[i], self.location.lon, u.deg, 11)
            self.assertQuantity(r.height[i], self.location.height, u.m, 6)

        # Check the vetcorized Cartesian conversion
        r = GeodeticRepresentation(latitude, longitude, height).to_cartesian()

        for i in range(n):
            self.assertQuantity(r.x[i], cart.x[i], u.m, 6)
            self.assertQuantity(r.y[i], cart.y[i], u.m, 6)
            self.assertQuantity(r.z[i], cart.z[i], u.m, 6)


    def test_horizontal(self):
        for (angle, point) in (((90, 0), (1, 0, 0)),
                               (( 0, 0), (0, 1, 0)),
                               (( 0, 90), (0, 0, 1)),
                               (( -90, 0), (-1, 0, 0))):
            h = HorizontalRepresentation(azimuth=angle[0] * u.deg,
                                         elevation=angle[1] * u.deg)
            cart = h.represent_as(CartesianRepresentation)

            self.assertQuantity(cart.x, point[0], u.one, 9)
            self.assertQuantity(cart.y, point[1], u.one, 9)
            self.assertQuantity(cart.z, point[2], u.one, 9)

            cart = CartesianRepresentation(*point)
            h = cart.represent_as(HorizontalRepresentation)

            self.assertQuantity(h.azimuth, angle[0] * u.deg, u.deg, 7)
            self.assertQuantity(h.elevation, angle[1] * u.deg, u.deg, 7)


    def test_ecef(self):
        # Check the forward transform
        ecef = ECEF(self.location.itrs.cartesian, obstime=self.obstime)
        itrs = ecef.transform_to(ITRS(obstime=self.obstime))
        self.assertCartesian(ecef, itrs, u.m, 8)

        # Check the backward transform
        ecef0 = itrs.transform_to(ECEF(obstime=self.obstime))
        self.assertCartesian(ecef, ecef0, u.m, 8)

        # Check the obstime handling
        ecef1 = itrs.transform_to(ECEF)
        self.assertEqual(ecef1.obstime, itrs.obstime)
        self.assertCartesian(ecef1, ecef0, u.m, 8)

        # Check the round trip with different obstimes
        itrs = ecef.transform_to(ITRS)
        ecef0 = itrs.transform_to(ECEF(obstime=self.obstime))
        self.assertCartesian(ecef, ecef0, u.m, 2)


    def test_enu(self):
        ecef = ECEF(self.location.itrs.cartesian, obstime=self.obstime)

        # Check the constructor & to ECEF transform
        enu = ENU(x=0 * u.m, y=0 * u.m, z=0 * u.m, location=self.location) 
        r = enu.transform_to(ECEF)

        self.assertEqual(r.obstime, enu.obstime)
        self.assertQuantity(r.x, ecef.x, u.m, 6)
        self.assertQuantity(r.y, ecef.y, u.m, 6)
        self.assertQuantity(r.z, ecef.z, u.m, 6)

        # Check the from ECEF transform
        enu = ecef.transform_to(ENU(location=self.location))

        self.assertEqual(enu.obstime, self.obstime)
        self.assertQuantity(enu.x, 0, u.m, 6)
        self.assertQuantity(enu.y, 0, u.m, 6)
        self.assertQuantity(enu.z, 0, u.m, 6)

        # Check the affine transform
        points = ((0, 0, 1), (1, 0, 0), (0, 1, 0), (1, 1, 0), (1, 1, 1),
                  (0, 1, 1))
        for point in (points):
            cart = CartesianRepresentation(x=point[1], y=point[0], z=point[2],
                                           unit=u.m)
            altaz = SkyCoord(cart, frame="altaz", location=self.location,
                             obstime=self.obstime)
            ecef0 = altaz.transform_to(ECEF(obstime=self.obstime))

            cart = CartesianRepresentation(x=point[0], y=point[1], z=point[2],
                                           unit=u.m)
            enu = ENU(cart, location=self.location, obstime=self.obstime)
            ecef1 = enu.transform_to(ECEF)

            self.assertEqual(ecef0.obstime, ecef1.obstime)
            self.assertQuantity(ecef0.x, ecef1.x, u.m, 4)
            self.assertQuantity(ecef0.y, ecef1.y, u.m, 4)
            self.assertQuantity(ecef0.z, ecef1.z, u.m, 4)

        # Check the orientation
        point = (1, -1, 2)
        cart = CartesianRepresentation(x=point[0], y=point[1], z=point[2],
                                       unit=u.m)
        altaz = SkyCoord(cart, frame="altaz", location=self.location,
                         obstime=self.obstime)
        ecef0 = altaz.transform_to(ECEF(obstime=self.obstime))

        for (orientation, sign) in ((("N", "E", "U"), (1, 1, 1)),
                                    (("N", "E", "D"), (1, 1, -1)),
                                    (("S", "E", "U"), (-1, 1, 1)),
                                    (("N", "W", "U"), (1, -1, 1))):
            cart = CartesianRepresentation(x=sign[0] * point[0],
                y=sign[1] * point[1], z=sign[2] * point[2], unit=u.m)
            enu = ENU(cart, location=self.location, obstime=self.obstime,
                      orientation=orientation)
            ecef1 = enu.transform_to(ECEF(obstime=self.obstime))

            self.assertQuantity(ecef0.x, ecef1.x, u.m, 4)
            self.assertQuantity(ecef0.y, ecef1.y, u.m, 4)
            self.assertQuantity(ecef0.z, ecef1.z, u.m, 4)

        # Check the unit vector case
        uy = HorizontalRepresentation(azimuth = 0 * u.deg,
                                      elevation = 0 * u.deg)
        enu = ENU(uy, location=self.location, obstime=self.obstime)

        self.assertQuantity(enu.x, 0, u.one, 9)
        self.assertQuantity(enu.y, 1, u.one, 9)
        self.assertQuantity(enu.z, 0, u.one, 9)

        r = enu.transform_to(ECEF)

        self.assertEqual(r.obstime, enu.obstime)
        self.assertQuantity(r.cartesian.norm(), 1, u.one, 6)

        ecef = ECEF(uy, obstime=self.obstime)
        enu = ecef.transform_to(ENU(location=self.location))

        self.assertEqual(enu.obstime, ecef.obstime)
        self.assertQuantity(enu.cartesian.norm(), 1, u.one, 6)

        # Check the magnetic north case
        enu0 = ENU(uy, location=self.location, obstime=self.obstime)
        frame1 = ENU(location=self.location, obstime=self.obstime,
                     magnetic=True)
        enu1 = enu0.transform_to(frame1)
        self.assertEqual(enu0.obstime, enu1.obstime)

        declination = numpy.arcsin(enu0.cartesian.cross(enu1.cartesian).norm())
        self.assertQuantity(declination.to(u.deg), 1.11 * u.deg, u.deg, 2)

        # Test the magnetic case with no obstime
        with self.assertRaises(ValueError) as context:
            ENU(uy, location=self.location, magnetic=True)
        self.assertRegex(context.exception.args[0], "^Magnetic")

        # Test the invalid frame case
        with self.assertRaises(ValueError) as context:
            ENU(uy, location=self.location, orientation=("T", "O", "T", "O"))
        self.assertRegex(context.exception.args[0], "^Invalid frame")

        # Test the enu round trip with a position
        frame0 = ENU(location=self.location)
        frame1 = ENU(location=self.location, obstime=self.obstime,
                     magnetic=True)
        enu0 = ENU(x=1 * u.m, y=2 * u.m, z=3 * u.m, location=self.location,
                   obstime=self.obstime)
        enu1 = enu0.transform_to(frame1).transform_to(frame0)
        self.assertCartesian(enu0, enu1, u.m, 8)

        # Test the same frame case
        enu1 = enu0.transform_to(frame0)
        self.assertCartesian(enu0, enu1, u.m, 8)
        self.assertEqual(enu0.obstime, enu1.obstime)


if __name__ == "__main__":
    unittest.main()
