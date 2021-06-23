import numpy as np

A = 17.269
B = 35.86
C = 610.78
D = 273.16
K = 17.67
KK = 243.5
A0 = 6.107799961
A1 = 4.436518521e-1
A2 = 1.428945805e-2
A3 = 2.650648471e-4
A4 = 3.031240396e-6
A5 = 2.034080948e-8
A6 = 6.136820929e-11


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


# goff-gratch formula for water saturation vapour pressure
def vapor(tfp):
    t = tfp
    e = (
        -7.90298 * (373.16 / t - 1.0)
        + 5.02808 * np.log10(373.16 / t)
        - 1.3816e-7 * (10.0 ** (11.344 * (1.0 - t / 373.16)) - 1.0)
        + 8.1328e-3 * (10.0 ** (3.49149 * (1.0 - 373.16 / t)) - 1.0)
    )
    es = 1013.246 * 10.0 ** e
    satr = 100.0 * es

    return satr


# teten's formula for water saturation vapor pressure
def teten(t):
    return C * np.exp(A * (t - D) / (t - B))


# Bolton's formula for water saturation (Bolton, Mon. Weather Rev, 1980)
def bolton(temperature):
    t = temperature - 273.15
    return 611.2 * np.exp(K * t / (t + KK))


# Lowe and Ficke's formula for saturation (see Pruppacher & Klett, p625)
def lf(temperature):
    t = temperature - 273.15
    return 100.0 * (A0 + t * (A1 + t * (A2 + t * (A3 + t * (A4 + t * (A5 + A6 * t))))))


qsat_methods = {
    "goff-gratch": vapor,
    "teten": teten,
    "bolton": bolton,
    "lowe-ficke": lf,
}


def saturation_vapor_pressure(temperature, method="goff-gratch"):
    return qsat_methods[method.lower()](temperature)
