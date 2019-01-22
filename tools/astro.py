# -*- coding: utf-8 -*-
"""
Common framework for GRAND packages

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

import sys

from astropy.coordinates import EarthLocation, AltAz, Angle, ICRS
from astropy.time import Time


class AstroConversion:
    """
    AstroConversion

    This class handles astronomical conversion, such as coordinate transformation from alt-az on local site to sky
    coordinates.

    @todo Convert from print to logger
    """
    @classmethod
    def __init__(cls, longitude=None, latitude=None, altitude=None, ra=None, dec=None):
        cls.longitude = longitude
        cls.latitude = latitude
        cls.altitude = altitude
        cls.ra = ra
        cls.dec = dec

        if cls.longitude is not None and cls.latitude is not None and cls.altitude is not None:
            cls.localsite = EarthLocation.from_geodetic(lon=cls.longitude, lat=cls.latitude, height=cls.altitude,
                                                        ellipsoid='WGS84')
        else:
            print('AstroConversion: can not create a local site')
            sys.exit(1)

    @classmethod
    def to_skycoord(cls, theta=None, phi=None, time=None, coordsys=ICRS):
        """
        Transforms an input direction given in a local site from alt-az system (theta, phi) to
        `~astropy.coordinates.SkyCoord` sky coordinates.

        GRAND convention: receiver convention: phi is oriented West of North, theta from zenith

             z=Up
             /\
             |
             |
             | theta
             |- /.
             | / .
             |/  .
             --------------> y=West
            / .  .
           / / . .
          /-     .
         /   phi
        |/
       x=North

        :param theta: local azimuth of the event
        :param phi: local altitude of the event
        :param time: time to use to compute the AltAz coordinates, instance of `~astropy.time.Time`
        :param coordsys: coordinate system to use, can be an instance of `~astropy.coordinates` such as ICRS,
                         FK5, etc...
        :return: sky coordinates in 'ICRS'
        """
        az = -Angle(phi)
        zenith = '90d'
        alt = Angle(zenith)-Angle(theta)
        time = Time(time)
        c = AltAz(az=az, alt=alt, obstime=time, location=cls.localsite)
        return c.transform_to(coordsys)

    @classmethod
    def localsiderealtime(cls, time):
        """
        Provide the local sidereal time
        @todo Implement computation of local sidereal time
        :param time: input time, in whatever time scale supported by `~astropy.time.Time`
        :return: local sidereal time
        """
        t = Time(time, location=cls.localsite)
        return t.sidereal_time('apparent')
