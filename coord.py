'''
# Coordiante conversions. -----------------------------------------------------
#
# HMS, DMS to decimal degrees.
# Decimal degrees to HMS, DMS.
#
# -----------------------------------------------------------------------------
'''

from __future__ import print_function, division

import numpy
import math


def hms_dms_dd(ra, dec, delimiter=" "):
    """Convert from HMS; DMS to DD."""

    try: 
        ra_dd, dec_dd = float(ra), float(dec)
    except ValueError:

        if ":" in ra:
            delimiter = ":"
        elif "h" in ra:
            ra  = ra.replace("h", " ").replace("m", " ").replace("s", " ")
            dec = dec.replace("d", " ").replace("m", " ").replace("s", " ")

        ra, dec = ra.split(delimiter), dec.split(delimiter)

        # RA:
        ra_hours_dd = float(ra[0]) * 15.
        ra_minutes_dd = float(ra[1]) * 15. / 60.
        ra_seconds_dd = float(ra[2]) * 15. / 3600.
        ra_dd = ra_hours_dd + ra_minutes_dd + ra_seconds_dd
        if ra_dd >= 360.:
            ra_dd = abs(ra_dd  - 360.)

        # DEC:
        if "-" in dec[0]:
            dec_dd = float(dec[0]) - (float(dec[1]) / 60.) - (float(dec[2]) / 3600.)
        else:
            dec_dd = float(dec[0]) + (float(dec[1]) / 60.) + (float(dec[2]) / 3600.)
    

    return ra_dd, dec_dd




def HMS_DMS_DD(ra, dec, delimiter=" "):
    '''Converts right ascension and declination from hh:mm:ss, dd:mm:ss to
    decimal degrees (dd).

    Parameters
    ----------
    ra  : string
        Right Ascension (HMS). Can be delimited with any single character.
    dec : string
        Delination (DMS). Can be delimited with any single character.

    Returns
    -------
    ra_dd  : float
            Right ascension in decimal degrees. No rounding.
    dec_dd : float
            Declination in decimal degrees. No rounding.

    Raises
    ------
    ValueError
        If input RA and DEC are not in the correct format.

    Examples
    --------
    >>> a, b = hms_dms_dd('23 41 29', '-29 19 15')
    >>> a
    355.37083333333334
    >>> b
    -29.3208333333
    >>> c, d = hms_dms_dd('23:41:29', '29:19:15')
    >>> c
    355.37083333333334
    >>> d
    29.3208333333

    '''

    if ":" in ra:
        delimiter = ":"

    try:

        # Read RA string and split into parts:
        ra_hour = ra[0]+ra[1]
        ra_minute = ra[3]+ra[4]
        ra_second = ''
        for i in range(6, len(ra)):
            ra_second += ra[i]

        ra_hour = float(ra_hour)
        ra_minute = float(ra_minute)
        ra_second = float(ra_second)

        # Read DEC string and split into parts:
        if '-' in dec:
            dec_sign = '-'
            i = 1
        elif '+' in dec:
            dec_sign = '+'
            i = 1
        else:
            dec_sign = '+'
            i = 0

        dec_deg = dec[i]+dec[i+1]
        dec_arcmin = dec[i+3]+dec[i+4]
        dec_arcsec = ''
        for j in range(i+6, len(dec)):
            dec_arcsec += dec[j]

        dec_deg = float(dec_deg)
        dec_arcmin = float(dec_arcmin)
        dec_arcsec = float(dec_arcsec)

        # First convert RA from HMS to DMS:
        ra_deg = 15 * ra_hour
        ra_arcmin = 15 * ra_minute
        ra_arcsec = 15 * ra_second

        # DMS to decimal:
        # RA first:
        ra_deg += math.floor(ra_arcmin/60)
        ra_arcmin = numpy.remainder(ra_arcmin, 60)
        ra_arcmin += math.floor(ra_arcsec/60)
        ra_arcsec = numpy.remainder(ra_arcsec, 60)
        if ra_arcmin > 60:
            print('OOPS')  # This shouldn't happen. If it does: oops.
        ra_arcmin += numpy.true_divide(ra_arcsec, 60)
        ra_deg += numpy.true_divide(ra_arcmin, 60)

        # DEC second:
        dec_arcmin += numpy.true_divide(dec_arcsec, 60)
        dec_deg += numpy.true_divide(dec_arcmin, 60)
        dec_deg = str(dec_sign)+str(dec_deg)

        ra_dd, dec_dd = float(ra_deg), float(dec_deg)

        return ra_dd, dec_dd

    except ValueError:

        print('ValueError: Check input coordinates are in the right format.' \
            '\n\nFormat must be:' \
            '\nRight ascension [HMS]: hh:mm:ss.s or hh mm ss.s' \
            '\nDeclination     [DMS]: dd:mm:ss.s or dd mm ss.s')
        return None, None


def dd_hms_dms(ra, dec, delimiter=' '):
    '''Converts from decimal degress to HMS, DMS coordinates.

    Parameters
    ----------
    ra        : float
              Right ascension (decimal degrees).
    dec       : float
              Delcination (decimal degrees).
    delimiter : string, optional
              Delimiter character between the hours-minutes-second and
              degrees-arcmins-arcsecs. Default is a blank space.
    '''

    flag = False

    ra_h1 = numpy.true_divide(ra, 15.0)
    ra_h2 = math.floor(ra_h1)
    ra_m1 = 60 * (ra_h1 - ra_h2)
    ra_m2 = math.floor(ra_m1)
    ra_s = round(60 * (ra_m1 - ra_m2), 2)

    dec_d1 = dec
    if -1.0 <= dec_d1 < 0.0:  # To make sure the sign is carried through.
        dec_d2 = math.ceil(dec_d1)
        dec_m1 = abs(60 * (dec_d1 - dec_d2))
        flag = True
    elif dec_d1 < -1.0:
        dec_d2 = math.ceil(dec_d1)
        dec_m1 = abs(60 * (dec_d1 - dec_d2))
    elif dec_d1 >= 0.0:
        dec_d2 = math.floor(dec_d1)
        dec_m1 = (60 * (dec_d1 - dec_d2))
    dec_m2 = math.floor(dec_m1)
    dec_s = round(60 * (dec_m1 - dec_m2), 2)

    ra_h = int(ra_h2)
    ra_m = int(ra_m2)

    dec_d = int(dec_d2)
    dec_m = int(dec_m2)

    if ra_h < 10:
        ra_h = str('0{0}'.format(ra_h))
    else:
        ra_h = str(ra_h)
    if ra_m < 10:
        ra_m = str('0{0}'.format(ra_m))
    else:
        ra_m = str(ra_m)
    if ra_s < 10:
        ra_s = str('0{0}'.format(ra_s))
    else:
        ra_s = str(ra_s)

    if dec_d < 0:
        if dec_d > -10:
            dec_d = '-0' + str(dec_d)[1]
        else:
            dec_d = str(dec_d)
    else:
        if dec_d < 10:
            if flag == True:
                dec_d = '-0' + str(dec_d)[0]
            else:
                dec_d = '+0' + str(dec_d)[0]
        else:
            dec_d = str(dec_d)
    if dec_m < 10:
        dec_m = str('0{0}'.format(dec_m))
    else:
        dec_m = str(dec_m)
    if dec_s < 10:
        dec_s = str('0{0}'.format(dec_s))
    else:
        dec_s = str(dec_s)

    ra_hms = ra_h + delimiter + ra_m + delimiter + ra_s
    dec_dms = dec_d + delimiter + dec_m + delimiter + dec_s

    return ra_hms, dec_dms
