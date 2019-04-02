# -*- coding: utf-8 -*-
"""
Unit tests for the grand_tools.astro module
"""

import unittest

from astropy.coordinates import EarthLocation
from astropy.time import Time

import grand_tools


class AstroTest(unittest.TestCase):
    """Unit tests for the astro module"""

    def test_to_skycoord(self):
        t = Time('2001-01-01 12:00:00', scale='utc')
        g = EarthLocation.of_site('greenwich')
        a = grand_tools.astro.AstroConversion(longitude=g.lon,
                                              latitude=g.lat,
                                              altitude=g.height)
        c = a.to_skycoord(theta='0d', phi='0d', time=t)
        self.assertAlmostEqual(c.ra.value, 281.20741528)
        self.assertAlmostEqual(c.dec.value, 51.47647952)

    def test_localsiderealtime(self):
        a = grand_tools.astro.AstroConversion(longitude='0d',
                                              latitude='0d',
                                              altitude='0m')
        t = Time('2000-01-01T12:00:00', format='isot', scale='tt')
        self.assertAlmostEqual(a.localsiderealtime(t).value, 18.67935935)


if __name__ == "__main__":
    unittest.main()
