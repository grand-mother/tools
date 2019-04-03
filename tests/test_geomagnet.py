# -*- coding: utf-8 -*-
"""
Unit tests for the grand_tools.geomagnet module
"""

import unittest

import grand_tools
from grand_tools.coordinates import ENU, GeodeticRepresentation
from grand_tools.geomagnet import Geomagnet

import numpy
import astropy.units as u
from astropy.coordinates import EarthLocation, ITRS


class GeomagnetTest(unittest.TestCase):
    """Unit tests for the geomagnet module"""


    def __init__(self, *args):
        super().__init__(*args)

        # The geo-magnetic field according to
        # http://geomag.nrcan.gc.ca/calc/mfcal-en.php"""
        self.ref = (0, 2.2983E-05, -4.0852E-05)

        # The corresponding Earth location
        self.location = EarthLocation(lat=45.0 * u.deg, lon=3.0 * u.deg,
                                      height=1000. * u.m)
        self.date = "2018-06-04"


    def assertField(self, field, tol=6):
        """Check that the magnetic field is consistent"""
        self.assertAlmostEqual((field.x / u.T).value, self.ref[0], tol)
        self.assertAlmostEqual((field.y / u.T).value, self.ref[1], tol)
        self.assertAlmostEqual((field.z / u.T).value, self.ref[2], tol)

    def get_coordinates(self, n=1):
        if n == 1:
            return ENU(x=0 * u.m, y=0 * u.m, z=0 * u.m, location=self.location,
                       obstime=self.date)
        else:
            zero = n * (0 * u.m,)
            return ENU(x=zero, y=zero, z=zero, location=self.location,
                       obstime=self.date)

    def test_default(self):
        # Test the initialisation
        model = "IGRF12"
        self.assertEqual(grand_tools.geomagnet.model(), model)
        self.assertEqual(grand_tools.geomagnet._default_magnet, None)

        # Test the default field getter
        c = self.get_coordinates()
        field = grand_tools.geomagnet.field(c)
        self.assertEqual(field.x.size, 1)
        self.assertEqual(c.obstime, field.obstime)
        self.assertField(field)

        # Test the vectorized getter
        n = 10
        c = self.get_coordinates(n)
        field = grand_tools.geomagnet.field(c)
        self.assertEqual(field.x.size, n)
        self.assertEqual(c.obstime, field.obstime)
        for value in field:
            frame = ENU(location=self.location)
            enu = value.transform_to(frame)
            self.assertField(enu)


    def test_custom(self):
        magnet = Geomagnet(model="WMM2015")
        c = self.get_coordinates()
        field = magnet(c)
        self.assertEqual(field.x.size, 1)
        self.assertEqual(c.obstime, field.obstime)
        self.assertField(field)


if __name__ == "__main__":
    unittest.main()
