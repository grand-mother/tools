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
    """Get the geoid undulation

    This method accepts vectorized input. The returned value is in unit
    `astropy.Unit.m`.

    Parameters
    ----------
    coordinates : ECEF or ENU
        The coordinates of points where the undulation is requested

    Returns
    -------
    astropy.Quantity
        The geoid undulation at the given point(s)


    Examples
    --------
    ```
    >>> from grand_tools.coordinates import ECEF
    >>> from grand_tools import topography

    >>> coordinates = ECEF(representation_type="geodetic", latitude=45,
    ...                    longitude=3, obstime="2019-01-01")
    >>> undulation = topography.geoid_undulation(coordinates)

    ```
    """

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

    def __init__(self, path):
        """Initialise a topography wrapper

        Alternative geomagnetic models can be instanciated with this class.
        Instead, use the `field` function of the geomagnet module in order to
        get the GRAND standard geomagnetic model.

        Parameters
        ----------
        path : str
            Path to the topography data

        Raises
        ------
        grand_libs.turtle.LibraryError
            The provided data are not valid / available

        Examples
        --------
        ```
        >>> from grand_tools.topography import Topography

        >>> topography = Topography("SRTMGL1")

        ```
        """
        self._stack = _Stack(path)


    def elevation(self, coordinates):
        """
        Get the topography elevation, w.r.t. sea level

        This method accepts vectorized input. The returned value is in unit
        `astropy.Unit.m`.

        Parameters
        ----------
        coordinates : ECEF or ENU
            The coordinates of points where the elevation is requested

        Returns
        -------
        astropy.Quantity
            The topography elevation at the given point(s)


        Examples
        --------
        ```
        >>> from grand_tools.coordinates import ECEF
        >>> from grand_tools.topography import Topography

        >>> topography = Topography("SRTMGL1")
        >>> coordinates = ECEF(representation_type="geodetic", latitude=45,
        ...                    longitude=3, obstime="2019-01-01")
        >>> elevation = topography.elevation(coordinates)

        ```
        """
        # Compute the geodetic coordinates
        cart = coordinates.transform_to(ECEF).cartesian
        geodetic = cart.represent_as(GeodeticRepresentation)

        # Return the topography elevation
        z = self._stack.elevation((geodetic.latitude / u.deg).value,
                                  (geodetic.longitude / u.deg).value)
        return z * u.m
