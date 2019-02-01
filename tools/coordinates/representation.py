# -*- coding: utf-8 -*-
"""
GRAND extension of astropy.coordinates

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

from collections import OrderedDict

from shared_libs import turtle

import numpy
import astropy.units as u
from astropy.coordinates import BaseRepresentation, CartesianRepresentation

__all__ = ["GeodeticRepresentation"]


class GeodeticRepresentation(BaseRepresentation):
    """Geodetic coordinates representation assuming WGS84 ellipsoid"""

    attr_classes = OrderedDict([["latitude", u.Quantity],
                                ["longitude", u.Quantity],
                                ["height", u.Quantity]])


    @classmethod
    def from_cartesian(cls, cart):
        m1 = 1 / u.m
        x, y, z = map(lambda v: v * m1, (cart.x, cart.y, cart.z))
        if x.size > 1:
            ecef = numpy.column_stack((x.value, y.value, z.value))
        else:
            ecef = (x.value, y.value, z.value)

        geodetic = turtle.ecef_to_geodetic(ecef)
        return cls(geodetic[0] * u.deg, geodetic[1] * u.deg, geodetic[2] * u.m)


    def to_cartesian(self):
        d1, m1 = 1 / u.deg, 1 / u.m
        ecef = turtle.ecef_from_geodetic(self.latitude * d1,
                                         self.longitude * d1, self.height * m1)
        if ecef.size == 3:
            return CartesianRepresentation(ecef[0] * u.m, ecef[1] * u.m,
                                           ecef[2] * u.m)
        else:
            return CartesianRepresentation(ecef[:,0] * u.m, ecef[:,1] * u.m,
                                           ecef[:,2] * u.m)
