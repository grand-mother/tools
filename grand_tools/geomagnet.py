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

from grand_libs.gull import Snapshot as _Snapshot
from grand_tools.coordinates import ECEF, ENU, GeodeticRepresentation

import numpy

import astropy.units as u
from astropy.coordinates import EarthLocation


_DEFAULT_MODEL = "IGRF12"
"""The default geo-magnetic model"""


_default_magnet = None
"""An instance of Geomagnet with the default geo-magnetic model"""


class Geomagnet:
    """Proxy to a geomagnetic model"""
    def __init__(self, model=None):
        if model is None:
            model = _DEFAULT_MODEL
        self._model = model
        self._snapshot = None
        self._date = None

    def __call__(self, coordinates):
        """
        Raises
        ------
        ValueError
            The provided coordinates are not valid
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
    """Get the geo-magnetic field"""
    global _default_magnet

    if _default_magnet is None:
        _default_magnet = Geomagnet()
    return _default_magnet(coordinates)


def model():
    """Get the default model for the geo-magnetic field"""
    return _DEFAULT_MODEL
