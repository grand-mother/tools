# -*- coding: utf-8 -*-
"""
Topography wrapper for GRAND packages

Copyright (C) 2018 The GRAND collaboration

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>
"""

from . import DATADIR
from .coordinates import ECEF, GeodeticRepresentation
from grand_libs.turtle import Stack as _Stack, Map as _Map

import os

import astropy.units as u
import numpy


_geoid = None
"""Map with geoid undulations"""


def geoid_undulation(coordinates):
    """Get the geoid undulation"""

    global _geoid

    if _geoid is None:
        path = os.path.join(DATADIR, "egm96.png")
        _geoid = _Map(path)

    # Compute the geodetic coordinates
    cart = coordinates.transform_to(ECEF).cartesian
    geodetic = cart.represent_as(GeodeticRepresentation)

    z = _geoid.elevation((geodetic.longitude / u.deg).value,
                         (geodetic.latitude / u.deg).value)
    return z * u.m


class Topography:
    """Proxy to topography data"""

    def __init__(self, path=None):
        if path is None:
            path = _DEFAULT_PATH
        self._stack = _Stack(path)

    def elevation(self, coordinates):
        # Compute the geodetic coordinates
        cart = coordinates.transform_to(ECEF).cartesian
        geodetic = cart.represent_as(GeodeticRepresentation)

        # Return the topography elevation
        z = self._stack.elevation((geodetic.longitude / u.deg).value,
                                  (geodetic.latitude / u.deg).value)
        return z * u.m
