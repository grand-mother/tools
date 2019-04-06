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

from .representation import GeodeticRepresentation, HorizontalRepresentation

from grand_libs import turtle

import numpy
import astropy.units as u
from astropy.coordinates import Attribute, BaseCoordinateFrame,                \
                                CartesianRepresentation,                       \
                                CylindricalRepresentation, FunctionTransform,  \
                                ITRS, PhysicsSphericalRepresentation,          \
                                RepresentationMapping, TimeAttribute,          \
                                frame_transform_graph


__all__ = ["ECEF", "ENU"]


_HAS_GEOMAGNET = False
"""The geomagnet module needs a deferred import due to circular references"""


class ECEF(BaseCoordinateFrame):
    """Earth-Centered Earth-Fixed frame, co-moving with the Earth
    """

    default_representation = CartesianRepresentation
    """Default representation of local frames"""

    obstime = TimeAttribute(default=None)
    """The observation time"""


    def __init__(self, *args, obstime=None, **kwargs):
        """Initialisation of an ECEF frame

        Parameters
        ----------
        obstime : Time or datetime or str, optional
            The observation time
        *args
            Any representation of the frame data, e.g. x, y, and z coordinates
        **kwargs
            Any extra BaseCoordinateFrame arguments
        """
        super().__init__(*args, obstime=obstime, **kwargs)


@frame_transform_graph.transform(FunctionTransform, ITRS, ECEF)
def itrs_to_ecef(itrs, ecef):
    """Compute the transformation from ITRS to ECEF coordinates

    Parameters
    ----------
    itrs : ITRS
        The initial coordinates in ITRS
    ecef : ECEF
        The ECEF frame to transform to

    Returns
    -------
    ECEF
        The ECEF frame with transformed coordinates
    """
    if ecef._obstime == None:
        # The conversion goes to a generic frame
        ecef._obstime = itrs._obstime
    elif ecef._obstime != itrs._obstime:
        itrs = itrs.transform_to(ITRS(obstime=ecef._obstime))

    return ecef.realize_frame(itrs.cartesian)


@frame_transform_graph.transform(FunctionTransform, ECEF, ITRS)
def ecef_to_itrs(ecef, itrs):
    """Compute the transformation from ECEF to ITRS coordinates

    Parameters
    ----------
    ecef : ECEF
        The initial coordinates in ECEF
    itrs : ITRS
        The ITRS frame to transform to

    Returns
    -------
    ITRS
        The ITRS frame with transformed coordinates
    """
    if itrs._obstime == ecef._obstime:
        c = ecef.cartesian
    else:
        itrs0 = ITRS(ecef.cartesian, obstime=ecef._obstime)
        c = itrs0.transform_to(ITRS(obstime=itrs._obstime)).cartesian

    return itrs.realize_frame(c)


@frame_transform_graph.transform(FunctionTransform, ECEF, ECEF)
def ecef_to_ecef(ecef0, ecef1):
    """Compute the transformation from ECEF to ECEF coordinates

    Parameters
    ----------
    ecef0 : ECEF
        The initial coordinates in the 1st ECEF frame
    ecef1 : ECEF
        The 2nd ECEF frame to transform to

    Returns
    -------
    ECEF
        The ECEF frame with transformed coordinates
    """
    if ecef1._obstime is None:
        ecef1._obstime = ecef0._obstime

    return ecef1.realize_frame(ecef0.cartesian)


class ENU(BaseCoordinateFrame):
    """Local geographic frames on the Earth, oriented along cardinal directions
    """

    default_representation = CartesianRepresentation
    """Default representation of local frames"""

    location = Attribute(default=None)
    """The origin on Earth of the local frame"""

    orientation = Attribute(default=("E", "N", "U"))
    """The orientation of the local frame, as cardinal directions"""

    magnetic = Attribute(default=False)
    """When enabled, use the magnetic North instead of the geographic one"""

    obstime = TimeAttribute(default=None)
    """The observation time"""


    def __init__(self, *args, location=None, orientation=None, magnetic=False,
                 obstime=None, **kwargs):
        """Initialisation of a local frame

        Parameters
        ----------
        *args
            Any representation of the frame data, e.g. x, y, and z coordinates
        location : EarthLocation
            The location on Earth of the local frame origin
        orientation : sequence of str, optional
            The cardinal directions of the x, y, and z axis (default: E, N, U)
        magnetic : boolean, optional
            Use the magnetic north instead of the geographic one (default: false)
        obstime : Time or datetime or str, optional
            The observation time
        **kwargs
            Any extra BaseCoordinateFrame arguments

        Raises
        ------
        ValueError
            The local frame configuration is not valid
        """

        # Do the base initialisation
        location = self.location if location is None else location
        orientation = self.orientation if orientation is None else orientation

        super().__init__(*args, location=location, orientation=orientation,
                         magnetic=magnetic, obstime=obstime, **kwargs)

        # Set the transform parameters
        itrs = self._location.itrs
        geo = itrs.represent_as(GeodeticRepresentation)
        latitude, longitude = geo.latitude / u.deg, geo.longitude / u.deg

        if magnetic:
            # Compute the magnetic declination
            if self._obstime is None:
                raise ValueError("Magnetic coordinates require specifying "
                                 "an observation time")
            ecef = ECEF(itrs.x, itrs.y, itrs.z, obstime=self._obstime)

            if not _HAS_GEOMAGNET:
                from ..geomagnet import field as _geomagnetic_field
            field = _geomagnetic_field(ecef)

            c = field.cartesian
            c /= c.norm()
            h = c.represent_as(HorizontalRepresentation)
            azimuth0 = (h.azimuth / u.deg).value
        else:
            azimuth0 = 0.

        def vector(name):
            tag = name[0].upper()
            if tag == "E":
                return turtle.ecef_from_horizontal(latitude, longitude,
                                                   90 + azimuth0, 0)
            elif tag == "W":
                return turtle.ecef_from_horizontal(latitude, longitude,
                                                   270 + azimuth0, 0)
            elif tag == "N":
                return turtle.ecef_from_horizontal(latitude, longitude,
                                                   azimuth0,  0)
            elif tag == "S":
                return turtle.ecef_from_horizontal(latitude, longitude,
                                                   180 + azimuth0,  0)
            elif tag == "U":
                return turtle.ecef_from_horizontal(latitude, longitude, 0, 90)
            elif tag == "D":
                return turtle.ecef_from_horizontal(latitude, longitude, 0, -90)
            else:
                raise ValueError(f"Invalid frame orientation `{name}`")

        ux = vector(self._orientation[0])
        uy = vector(self._orientation[1])
        uz = vector(self._orientation[2])

        self._basis = numpy.column_stack((ux, uy, uz))
        self._origin = itrs.cartesian


@frame_transform_graph.transform(FunctionTransform, ECEF, ENU)
def itrs_to_enu(ecef, enu):
    """Compute the transformation from ECEF to ENU coordinates

    Parameters
    ----------
    ecef : ECEF
        The initial coordinates in ECEF
    enu : ENU
        The ENU frame to transform to

    Returns
    -------
    ENU
        The ENU frame with transformed coordinates
    """
    c = ecef.cartesian
    if c.x.unit.is_equivalent("m"):
        c = c.copy()
        c -= enu._origin
    c = c.transform(enu._basis.T)

    if enu._obstime is None:
        enu._obstime = ecef._obstime

    return enu.realize_frame(c)


@frame_transform_graph.transform(FunctionTransform, ENU, ECEF)
def enu_to_itrs(enu, ecef):
    """Compute the transformation from ENU to ECEF coordinates

    Parameters
    ----------
    enu : ENU
        The initial coordinates in ENU
    ecef : ECEF
        The ECEF frame to transform to

    Returns
    -------
    ECEF
        The ECEF frame with transformed coordinates
    """
    c = enu.cartesian.transform(enu._basis)
    if c.x.unit.is_equivalent("m"):
        c += enu._origin

    if ecef._obstime is None:
        ecef._obstime = enu._obstime

    return ecef.realize_frame(c)


@frame_transform_graph.transform(FunctionTransform, ENU, ENU)
def enu_to_enu(enu0, enu1):
    """Compute the transformation from ENU to ENU coordinates

    Parameters
    ----------
    enu0 : ENU
        The initial coordinates in the 1st ENU frame
    enu1 : ENU
        The 2nd ENU frame to transform to

    Returns
    -------
    ENU
        The ENU frame with transformed coordinates
    """
    c = enu0.cartesian
    translate = c.x.unit.is_equivalent("m")

    # Forward the observation time
    if enu1._obstime is None:
        enu1._obstime = enu0._obstime

    # Check if the two frames are identicals
    if numpy.array_equal(enu0._basis, enu1._basis):
        if not translate or ((enu0._origin.x == enu1._origin.x) and
                             (enu0._origin.y == enu1._origin.y) and
                             (enu0._origin.z == enu1._origin.z)):
            # CartesianRepresentations might not eveluate to equal though the
            # coordinates are equal
            return enu1.realize_frame(c)

    # Transform to ECEF
    if translate:
        c = c.copy()
        c -= enu0._origin
    c = c.transform(enu0._basis.T)

    # Transform back from ECEF
    c = c.transform(enu1._basis)
    if translate:
        c += enu1._origin

    return enu1.realize_frame(c)
