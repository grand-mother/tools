# -*- coding: utf-8 -*-
"""
Extra frame(s) for astropy.coordinates

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

from grand_libs import turtle

import numpy
import astropy.units as u
from astropy.coordinates import Attribute, BaseCoordinateFrame,                \
                                CartesianRepresentation,                       \
                                CylindricalRepresentation, FunctionTransform,  \
                                ITRS, PhysicsSphericalRepresentation,          \
                                RepresentationMapping, TimeAttribute,          \
                                frame_transform_graph
from astropy.coordinates.builtin_frames.utils import DEFAULT_OBSTIME



__all__ = ["ENU"]


class ENU(BaseCoordinateFrame):
    """Local geographic frames on the Earth, oriented along cardinal directions
    """

    default_representation = CartesianRepresentation
    """Default representation of local frames"""

    location = Attribute(default=None)
    """The origin on Earth of the local frame"""

    orientation = Attribute(default=("E", "N", "U"))
    """The orientation of the local frame, as cardinal directions"""

    obstime = TimeAttribute(default=DEFAULT_OBSTIME)
    """The observation time"""


    def __init__(self, *args, location=None, orientation=None, **kwargs):
        """Initialisation of a local frame

        Parameters
        ----------
        *args
            Any representation of the frame data, e.g. x, y, and z coordinates
        location : EarthLocation
            The location on Earth of the local frame origin
        orientation : sequence of str, optional
            The cardinal directions of the x, y, and z axis (default: E, N, U)
        **kwargs
            Any extra BaseCoordinateFrame arguments

        Raises
        ------
        ValueError
            The local frame orientation is not valid
        """

        # Do the base initialisation
        location = self.location if location is None else location
        orientation = self.orientation if orientation is None else orientation

        super().__init__(*args, location=location, orientation=orientation,
                         **kwargs)

        # Set the transform parameters
        itrs = self._location.itrs
        geo = itrs.represent_as(GeodeticRepresentation)
        latitude, longitude = geo.latitude / u.deg, geo.longitude / u.deg

        def vector(name):
            tag = name[0].upper()
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
            else:
                raise ValueError(f"invalid frame orientation `{name}`")

        ux = vector(self._orientation[0])
        uy = vector(self._orientation[1])
        uz = vector(self._orientation[2])

        self._basis = numpy.column_stack((ux, uy, uz))
        self._origin = itrs.cartesian


@frame_transform_graph.transform(FunctionTransform, ITRS, ENU)
def itrs_to_enu(itrs, enu):
    """Compute the transformation from ITRS to ENU coordinates

    Parameters
    ----------
    itrs : ITRS
        The initial coordinates in ITRS
    enu : ENU
        The ENU frame to transform to

    Returns
    -------
    ENU
        The ENU frame with transformed coordinates
    """
    c = itrs.cartesian
    if c.x.unit.is_equivalent("m"):
        c._x -= enu._origin.x
        c._y -= enu._origin.y
        c._z -= enu._origin.z
    r = c.transform(enu._basis.T)

    return enu.realize_frame(r)


@frame_transform_graph.transform(FunctionTransform, ENU, ITRS)
def enu_to_itrs(enu, itrs):
    """Compute the transformation from ENU to ITRS coordinates

    Parameters
    ----------
    enu : ENU
        The initial coordinates in ENU
    itrs : ITRS
        The ITRS frame to transform to

    Returns
    -------
    ITRS
        The ITRS frame with transformed coordinates
    """
    r = enu.cartesian.transform(enu._basis)
    if r.x.unit.is_equivalent("m"):
        r._x += enu._origin.x
        r._y += enu._origin.y
        r._z += enu._origin.z
    r = CartesianRepresentation(r.x, r.y, r.z)

    return itrs.realize_frame(r)
