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

from .representation import GeodeticRepresentation 

from shared_libs import turtle

import numpy
import astropy.units as u
from astropy.coordinates import Attribute, BaseCoordinateFrame,                \
                                CartesianRepresentation,                       \
                                CylindricalRepresentation, FunctionTransform,  \
                                ITRS, PhysicsSphericalRepresentation,          \
                                RepresentationMapping, TimeAttribute,          \
                                frame_transform_graph

__all__ = ["ENU"]


class ENU(BaseCoordinateFrame):
    default_representation = CartesianRepresentation

    location = Attribute(default=None)
    """The origin of the local frame"""

    orientation = Attribute(default=("E", "N", "U"))
    """The orientation of the local frame"""

    equinox = TimeAttribute(default="B1950")
    obstime = TimeAttribute(default=None, secondary_attribute="equinox")


    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        # Set the transform parameters
        itrs = self._location.itrs
        geo = itrs.represent_as(GeodeticRepresentation)
        latitude, longitude = geo.latitude / u.deg, geo.longitude / u.deg

        def vector(tag):
            tag = tag[0].upper()
            if tag == "E":
                return turtle.ecef_from_horizontal(latitude, longitude, 90, 0)
            elif tag == "W":
                return turtle.ecef_from_horizontal(latitude, longitude, 270, 0)
            elif tag == "N":
                return turtle.ecef_from_horizontal(latitude, longitude, 0,  0)
            elif tag == "S":
                return turtle.ecef_from_horizontal(latitude, longitude, 180,  0)
            elif tag == "U":
                return turtle.ecef_from_horizontal(latitude, longitude, 0, 90)
            elif tag == "D":
                return turtle.ecef_from_horizontal(latitude, longitude, 0, -90)

        ux = vector(self._orientation[0])
        uy = vector(self._orientation[1])
        uz = vector(self._orientation[2])

        self._basis = numpy.column_stack((ux, uy, uz))
        self._origin = itrs.cartesian


@frame_transform_graph.transform(FunctionTransform, ITRS, ENU)
def itrs_to_enu(itrs, enu):
    """Compute the transformation from ITRS to ENU coordinates"""
    c = itrs.cartesian
    c._x -= enu._origin.x
    c._y -= enu._origin.y
    c._z -= enu._origin.z
    r = c.transform(enu._basis.T)

    return enu.realize_frame(r)


@frame_transform_graph.transform(FunctionTransform, ENU, ITRS)
def enu_to_itrs(enu, itrs):
    """Compute the transformation from ENU to ITRS coordinates"""
    r = enu.cartesian.transform(enu._basis)
    r._x += enu._origin.x
    r._y += enu._origin.y
    r._z += enu._origin.z
    r = CartesianRepresentation(r.x, r.y, r.z)

    return itrs.realize_frame(r)
