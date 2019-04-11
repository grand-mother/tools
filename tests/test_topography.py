# -*- coding: utf-8 -*-
"""
Unit tests for the grand_tools.topography module
"""

import os
import unittest

import grand_store
import grand_tools
from grand_tools.coordinates import ECEF, GeodeticRepresentation
from grand_tools.topography import Topography

import numpy
import astropy.units as u


class TopographyTest(unittest.TestCase):
    """Unit tests for the topography module"""


    def test_geoid(self):
        # Test the undulation getter
        c = ECEF(GeodeticRepresentation(latitude=45.5, longitude=3.5))
        z = grand_tools.topography.geoid_undulation(c)
        self.assertEqual(z.size, 1)
        self.assertEqual(z.unit, u.m)


    def test_topography(self):
        # Fetch a test tile
        dirname, basename = "tests/topography", "N39E092.SRTMGL1.hgt"
        path = os.path.join(dirname, basename)
        if not os.path.exists(path):
            try:
                os.makedirs(dirname)
            except OSError:
                pass
            with open(path, "wb") as f:
                f.write(grand_store.get(basename))

        # Test the topography getter
        topo = Topography(dirname)
        c = ECEF(GeodeticRepresentation(latitude=39.5, longitude=92.5))
        z = topo.elevation(c)
        self.assertEqual(z.size, 1)
        self.assertEqual(z.unit, u.m)
        self.assertFalse(numpy.isnan(z))


if __name__ == "__main__":
    unittest.main()
