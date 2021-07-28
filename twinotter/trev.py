import numpy as np
from .derive import saturation_vapor_pressure


def trev(cloud_base_pressure, temperature, pressure, qsat_method="goff-gratch"):
    """
    calculates temperature t and liquid water content alwc from cloud base pressure po
    and temperature to, for adiabatic ascent to the pressure p

    pressure in Pa, T in K
    returns

    Args:
        cloud base pressure po (Pa) and temperature to (K) and pressure at observation
        level (Pa)

    Returns:
        "adiabatic" temperature (C) and liquid water content (g/m3)
    """
    eps = 0.622
    cpd = 1.0057e3
    cw = 4.18e3
    cpv = 1.875e3
    rd = 287.04
    alref = 2.43e06
    tref = 303.15

    # Saturation vapour pressure
    e = saturation_vapor_pressure(temperature, method=qsat_method)

    # Saturation mixing ratio
    r = eps * e / (cloud_base_pressure - e)

    # Total heat capacity?
    cpt = cpd + r * cw

    # alhv = alref + (cpv - cw) * (temperature - tref)
    alhv = (2.501 - 0.00237 * (temperature - 273.15)) * 1.0e6
    thetaq = (
        temperature
        * (1.0e5 / (cloud_base_pressure - e)) ** (rd / cpt)
        * np.exp(alhv * r / (cpt * temperature))
    )

    # 1st approx
    t1 = temperature
    e = saturation_vapor_pressure(t1, method=qsat_method)
    rv = eps * e / (pressure - e)
    t1 = thetaq / (
        (1.0e5 / (pressure - e)) ** (rd / cpt) * np.exp(alhv * rv / (cpt * t1))
    )

    # successive approximations
    for i in range(100):
        e = saturation_vapor_pressure(t1, method=qsat_method)
        rv = eps * e / (pressure - e)
        # alhv = alref + (cpv - cw) * (t1 - tref)
        alhv = (2.501 - 0.00237 * (t1 - 273.15)) * 1.0e6
        t1 = (
            thetaq
            / ((1.0e5 / (pressure - e)) ** (rd / cpt) * np.exp(alhv * rv / (cpt * t1)))
            + t1
        ) / 2.0

    t = t1 - 273.15

    # get lwc
    e = saturation_vapor_pressure(t1, method=qsat_method)
    rv = eps * e / (pressure - e)
    tw = r - rv
    alwc = tw * pressure * 28.9644 / (8.314e7 * t1) * 1.0e7

    return t, alwc, thetaq
