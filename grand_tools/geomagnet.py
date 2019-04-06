# -*- coding: utf-8 -*-
"""
Geomagnetic field wrapper for GRAND packages

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

from .coordinates import ECEF, ENU, GeodeticRepresentation
from grand_libs.gull import Snapshot as _Snapshot

import astropy.units as u
import numpy

from astropy.coordinates import EarthLocation


_DEFAULT_MODEL = "IGRF12"
"""The default geo-magnetic model"""


_default_magnet = None
"""An instance of Geomagnet with the default geo-magnetic model"""


class Geomagnet:
    """Proxy to a geomagnetic model"""

    def __init__(self, model=None):
        """Initialise a geomagnet wrapper

        Alternative geomagnetic models can be instanciated with this class.
        Instead, use the `field` function of the geomagnet module in order to
        get the GRAND standard geomagnetic model.

        Parameters
        ----------
        model : str, optional
            The geomagnetic model (defaults to "IGRF12")

        Raises
        ------
        grand_libs.gull.LibraryError
            The requested model is not valid / available

        Notes
        -----
        In order to get the geomagnetic field from the default GRAND model
        (IGRF12), one should **not** directly instantiate this class.  Instead
        one should use the `field` function of the geomagnet module.  This class
        is only meant to be used for studies of the impact of alternative
        models.

        Supported models are:
          1. [IGRF12](https://www.ngdc.noaa.gov/IAGA/vmod/igrf.html)
          2. [WMM2015](http://www.geomag.bgs.ac.uk/research/modelling/WorldMagneticModel.html)

        Examples
        --------
        ```
        >>> from grand_tools.geomagnet import Geomagnet

        >>> geomagnet = Geomagnet("WWM2015")

        ```
        """
        if model is None:
            model = _DEFAULT_MODEL
        self._model = model
        self._snapshot = None
        self._date = None

    def field(self, coordinates):
        """
        Get the geo-magnetic field components

        This method accepts vectorized input. The frame of the returned value
        depends on the input size. If a single point is provided the magnetic
        field is returned in local ENU coordinates, centered on the point.
        Otherwise the components are given in ECEF.

        Parameters
        ----------
        coordinates : ECEF or ENU
            The coordinates of points where the magnetic field is requested

        Returns
        -------
        ECEF or ENU
            The magnetic field components at the given point(s)

        Raises
        ------
        ValueError
            The provided coordinates are not valid

        Examples
        --------
        ```
        >>> from grand_tools.coordinates import ECEF
        >>> from grand_tools.geomagnet import Geomagnet

        >>> geomagnet = Geomagnet()
        >>> coordinates = ECEF(representation_type="geodetic", latitude=45,
        ...                    longitude=3, obstime="2019-01-01")
        >>> field = geomagnet.field(coordinates)

        ```
        """

        # Update the snapshot, if needed
        obstime = coordinates.obstime
        if obstime is None:
            raise ValueError(
                "No observation time was specified for the coordinates")
        date = obstime.datetime.date()
        if date != self._date:
            self._snapshot = _Snapshot(self._model, date)
            self._date = date

        # Compute the geodetic coordinates
        cart = coordinates.transform_to(ECEF).cartesian
        geodetic = cart.represent_as(GeodeticRepresentation)

        # Fetch the magnetic field components in local ENU
        field = self._snapshot(geodetic.latitude / u.deg,
                               geodetic.longitude / u.deg,
                               geodetic.height / u.m)

        # Encapsulate the result
        n = geodetic.latitude.size
        if n == 1:
            location = EarthLocation(lat=geodetic.latitude,
                                     lon=geodetic.longitude,
                                     height=geodetic.height)
            return ENU(x=field[0] * u.T, y=field[1] * u.T, z=field[2] * u.T,
                       location=location, obstime=obstime)
        else:
            ecef = numpy.zeros((n, 3))
            for i, value in enumerate(field):
                location = EarthLocation(lat=geodetic.latitude[i],
                                         lon=geodetic.longitude[i],
                                         height=geodetic.height[i])
                enu = ENU(x=value[0], y=value[1], z=value[2], location=location)
                c = enu.transform_to(ECEF).cartesian
                ecef[i,0] = c.x
                ecef[i,1] = c.y
                ecef[i,2] = c.z
            return ECEF(x=ecef[:,0] * u.T, y=ecef[:,1] * u.T, z=ecef[:,2] * u.T,
                        obstime=obstime)


def field(coordinates):
    """Get the geo-magnetic field

    This method accepts vectorized input. The frame of the returned value
    depends on the input size. If a single point is provided the magnetic
    field is returned in local ENU coordinates, centered on the point.
    Otherwise the components are given in ECEF.

    Parameters
    ----------
    coordinates : ECEF or ENU
        The coordinates of points where the magnetic field is requested

    Returns
    -------
    ECEF or ENU
        The magnetic field components at the given point(s)

    Raises
    ------
    ValueError
        The provided coordinates are not valid

    Examples
    --------
    ```
    >>> from grand_tools.coordinates import ECEF
    >>> from grand_tools import geomagnet

    >>> coordinates = ECEF(representation_type="geodetic", latitude=45,
    ...                    longitude=3, obstime="2019-01-01")
    >>> field = geomagnet.field(coordinates)

    ```
    """
    global _default_magnet

    if _default_magnet is None:
        _default_magnet = Geomagnet()
    return _default_magnet.field(coordinates)


def model():
    """Get the default model for the geo-magnetic field

    Currently the default geomagnetic model used in GRAND is
    [IGRF12](https://www.ngdc.noaa.gov/IAGA/vmod/igrf.html).

    Returns
    -------
    str
        The name of the default model

    Examples
    --------
    ```
    >>> from grand_tools import geomagnet

    >>> geomagnet.model()
    'IGRF12'

    ```
    """
    return _DEFAULT_MODEL
