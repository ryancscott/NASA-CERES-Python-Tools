
# This code is from the PyAstronomy package

import datetime
import numpy
import six.moves as smo


def daycnv(xjd, mode="dtlist"):
    """
    Converts Julian dates to Gregorian calendar dates.
    Handles both individual floats as xjd and iterables such as
    lists and arrays. In the latter case, the result is returned
    in the form of a list.
    Parameters
    ----------
    xjd : float, list, array
        The Julian date
    mode : string, {idl, dtlist, dt}, optional
        Determines format of output. If 'idl' is given (default),
        a list holding [year, month, day, (fractional) hours] is
        returned; this mimics the behavior of the IDL astrolib function.
        If 'dtlist' is given, a list holding
        [year, month, day, hours, minutes, seconds, microseconds] is
        returned. Finally, if 'dt' is specified, a Python
        datetime object will be returned. If the input is an iterable,
        the mode determines the format of the individual items in the
        result list.
    Returns
    -------
    Calendar date : list or datetime object
        A list holding [year, month, day, (fractional) hours] (default)
        or [year, month, day, hours, minutes, seconds, microseconds].
        Alternatively, a Python datetime object is returned. The format
        depends on the 'mode' specified. If the input is an iterable of
        Julian dates, the output is a list.
    Notes
    -----
    .. note:: This function was ported from the IDL Astronomy User's Library.
    :IDL - Documentation:
    NAME:
          DAYCNV
    PURPOSE:
          Converts Julian dates to Gregorian calendar dates
    CALLING SEQUENCE:
          DAYCNV, XJD, YR, MN, DAY, HR
    INPUTS:
          XJD = Julian date, positive double precision scalar or vector
    OUTPUTS:
          YR = Year (Integer)
          MN = Month (Integer)
          DAY = Day (Integer)
          HR = Hours and fractional hours (Real).   If XJD is a vector,
                  then YR,MN,DAY and HR will be vectors of the same length.
    EXAMPLE:
          IDL> DAYCNV, 2440000.D, yr, mn, day, hr
          yields yr = 1968, mn =5, day = 23, hr =12.
    WARNING:
          Be sure that the Julian date is specified as double precision to
          maintain accuracy at the fractional hour level.
    METHOD:
          Uses the algorithm of Fliegel and Van Flandern (1968) as reported in
          the "Explanatory Supplement to the Astronomical Almanac" (1992), p. 604
          Works for all Gregorian calendar dates with XJD > 0, i.e., dates after
          -4713 November 23.
    REVISION HISTORY:
          Converted to IDL from Yeoman's Comet Ephemeris Generator,
          B. Pfarr, STX, 6/16/88
          Converted to IDL V5.0   W. Landsman   September 1997
    """

    # if not mode in ('idl', 'dtlist', 'dt'):
    #     raise(PE.PyAValError("Unknown mode: " + str(mode),
    #                          where="daycnv",
    #                          solution="Use any of 'idl', 'dtlist', or 'dt'."))

    # Adjustment needed because Julian day starts at noon, calendar day at midnight

    iterable = hasattr(xjd, "__iter__")

    # Use iterable throughout calculations
    if not iterable:
        xjd = [xjd]

    jd = numpy.array(xjd).astype(int)  # Truncate to integral day
    frac = numpy.array(xjd).astype(float) - jd + \
        0.5  # Fractional part of calendar day
    gi = numpy.where(frac >= 1.0)
    frac[gi] -= 1.0
    jd[gi] += 1

    hr = frac * 24.0
    l = jd + 68569
    n = 4 * l // 146097
    l = l - (146097 * n + 3) // 4
    yr = 4000 * (l + 1) // 1461001
    l = l - 1461 * yr // 4 + 31  # 1461 = 365.25 * 4
    mn = 80 * l // 2447
    day = l - 2447 * mn // 80
    l = mn // 11
    mn = mn + 2 - 12 * l
    yr = 100 * (n - 49) + yr + l
    if mode in ('dt', 'dtlist'):
        # [year, month, day, hours, minutes, seconds, microseconds] requested
        hour = numpy.floor(hr).astype(int)
        minute = numpy.floor((hr - numpy.floor(hr)) * 60).astype(int)
        sec = numpy.floor((hr - hour - minute / 60.) * 3600.).astype(int)
        msec = (3600 * 1e6 * (hr - hour - minute / 60. - sec / 3600.)).astype(int)
        if mode == 'dtlist':
            if not iterable:
                return [yr[0], mn[0], day[0], hour[0], minute[0], sec[0], msec[0]]
            return [[yr[i], mn[i], day[i], hour[i], minute[i], sec[i], msec[i]] for i in smo.range(len(yr))]
        # Return datetime object
        dts = [datetime.datetime(*(yr[i], mn[i], day[i], hour[i],
                                   minute[i], sec[i], msec[i])) for i in smo.range(len(yr))]
        if not iterable:
            return dts[0]
        return dts
    if not iterable:
        return [yr[0], mn[0], day[0], hr[0]]
    return [[yr[i], mn[i], day[i], hr[i]] for i in smo.range(len(yr))]