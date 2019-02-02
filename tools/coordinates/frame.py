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
    CartesianRepresentation, FunctionTransform, ITRS, RepresentationMapping,   \
    TimeAttribute, frame_transform_graph

__all__ = ["ENU"]


class ENU(BaseCoordinateFrame):
    default_representation = CartesianRepresentation

    location = Attribute(default=None)
    """The origin of the ENU frame"""

    equinox = TimeAttribute(default="B1950")
    obstime = TimeAttribute(default=None, secondary_attribute="equinox")


    frame_specific_representation_info = {
        CartesianRepresentation: [
            RepresentationMapping("x", "x"),
            RepresentationMapping("y", "y"),
            RepresentationMapping("z", "z")]}


    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        # Set the transform parameters
        itrs = self._location.itrs
        geo = itrs.represent_as(GeodeticRepresentation)
        latitude, longitude = geo.latitude / u.deg, geo.longitude / u.deg
        ux = turtle.ecef_from_horizontal(latitude, longitude, 90, 0)
        uy = turtle.ecef_from_horizontal(latitude, longitude, 0, 90)
        uz = turtle.ecef_from_horizontal(latitude, longitude, 0, 90)

        self._basis = numpy.column_stack((ux, uy, uz))
        self._origin = itrs.cartesian


@frame_transform_graph.transform(FunctionTransform, ITRS, ENU)
def itrs_to_enu(itrs, enu):
    """Compute the transformation from ITRS to ENU coordinates"""
    c = itrs.cartesian
    c.x -= enu._origin.x
    c.y -= enu._origin.y
    c.z -= enu._origin.z
    r = c.transform(self._basis)

    return enu.realize_frame(r)


@frame_transform_graph.transform(FunctionTransform, ENU, ITRS)
def enu_to_itrs(enu, itrs):
    """Compute the transformation from ENU to ITRS coordinates"""
    r = enu.cartesian.transform(enu._basis.T)
    r = CartesianRepresentation(x = r.x + enu._origin.x,
        y = r.y + enu._origin.y, z = r.z + enu._origin.z)

    return itrs.realize_frame(r)
