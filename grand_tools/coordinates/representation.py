# -*- coding: utf-8 -*-
"""
Extra representations for astropy.coordinates

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

from grand_libs import turtle

import numpy
import astropy.units as u
from astropy.coordinates import BaseRepresentation, CartesianRepresentation

__all__ = ["GeodeticRepresentation", "HorizontalRepresentation"]


class GeodeticRepresentation(BaseRepresentation):
    """Geodetic coordinates representation w.r.t. the WGS84 ellipsoid"""

    attr_classes = OrderedDict([["latitude", u.Quantity],
                                ["longitude", u.Quantity],
                                ["height", u.Quantity]])
    """Attributes of a Geodetic representation"""


    def __init__(self, latitude, longitude, height, copy=True):
        """Initialise a geodetic representation

        Parameters
        ----------
        latitude : Quantity or str
            The latitude angle measured clockwise, w.r.t. the xOy plane
        longitude : Quantity or str
            The longitude angle measured counter-clockwise, w.r.t. the x-axis
        height : Quantity or str
            The height above the WGS84 ellipsoid

        copy : bool, optional
            If `True` (default), arrays will be copied rather than referenced
        """
        super().__init__(latitude, longitude, height, copy=copy)


    @classmethod
    def from_cartesian(cls, cart):
        """Generate a Geodetic representation from a Cartesian one

        Parameters
        ----------
        cart : CartesianRepresentation
            The cartesian coordinates of a point, e.g. in ITRS

        Returns
        -------
        GeodeticRepresentation
            The corresponding geodetic coordinates
        """
        m1 = 1 / u.m
        x, y, z = map(lambda v: v * m1, (cart.x, cart.y, cart.z))
        if x.size > 1:
            ecef = numpy.column_stack((x.value, y.value, z.value))
        else:
            ecef = (x.value, y.value, z.value)

        geodetic = turtle.ecef_to_geodetic(ecef)
        return cls(geodetic[0] * u.deg, geodetic[1] * u.deg, geodetic[2] * u.m,
                   copy=False)


    def to_cartesian(self):
        """Generate a Cartesian representation from a Geodetic one

        Returns
        -------
        CartesianRepresentation
            The Cartesian coordinates corresponding to this representation
        """
        d1, m1 = 1 / u.deg, 1 / u.m
        ecef = turtle.ecef_from_geodetic(self.latitude * d1,
                                         self.longitude * d1, self.height * m1)
        if ecef.size == 3:
            return CartesianRepresentation(ecef[0] * u.m, ecef[1] * u.m,
                                           ecef[2] * u.m, copy=False)
        else:
            return CartesianRepresentation(ecef[:,0] * u.m, ecef[:,1] * u.m,
                                           ecef[:,2] * u.m, copy=False)


class HorizontalRepresentation(BaseRepresentation):
    """Horizontal angular representation, for unit vectors"""


    attr_classes = OrderedDict([["azimuth", u.Quantity],
                                ["elevation", u.Quantity]])
    """Attributes of a Horizontal representation"""


    def __init__(self, azimuth, elevation, copy=True):
        """Initialise a Horizontal angular representation of a unit vector

        Parameters
        ----------
        azimuth : Quantity or str
            The azimuth angle measured clockwise, w.r.t. the y axis
        elevation : Quantity or str
            The elevation angle measured clockwise, w.r.t. the xOy plane

        copy : bool, optional
            If `True` (default), arrays will be copied rather than referenced
        """
        super().__init__(azimuth, elevation, copy=copy)


    @classmethod
    def from_cartesian(cls, cart):
        """Generate a Horizontal angular representation from a Cartesian unit
        vector

        **Note** that the Cartesian unit vector **must** be dimensioneless.
        Though it is not checked if the norm of the vector is indeed unity.

        Parameters
        ----------
        cart : CartesianRepresentation
            The cartesian coordinates of the unit vector

        Returns
        -------
        HorizontalRepresentation
            The corresponding angular coordinates

        Raises
        ------
        ValueError
            The cartesian representation is not dimensioneless
        """
        if ((cart.x.unit != u.one) or (cart.y.unit != u.one) or
            (cart.z.unit != u.one)):
            raise ValueError("coordinates must be dimensionless")

        rho = numpy.sqrt(cart.x**2 + cart.y**2)
        theta = numpy.arctan2(rho, cart.z)

        if theta == 0 * u.rad:
            elevation = 90 * u.deg
            azimuth = 0 * u.deg
        else:
            elevation = 90 * u.deg - theta
            azimuth = 90 * u.deg - numpy.arctan2(cart.y, cart.x)

        return cls(azimuth, elevation, copy=False)


    def to_cartesian(self):
        """Generate a Cartesian unit vector from this Horizontal angular
        representation

        Returns
        -------
        CartesianRepresentation
            The corresponding cartesian unit vector
        """
        theta = 90 * u.deg - self.elevation
        phi = 90 * u.deg - self.azimuth
        ct, st = numpy.cos(theta), numpy.sin(theta)

        return CartesianRepresentation(numpy.cos(phi) * st,
                                       numpy.sin(phi) * st, ct, copy=False)
