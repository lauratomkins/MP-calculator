import numpy as np
import constants


def eswFromTemp(temp_C):
    esw = 6.1094 * np.exp((17.625 * temp_C) / (243.04 + temp_C))

    return esw


def esiFromTemp(temp_C):
    esi = 6.1121 * np.exp((22.587 * temp_C) / (273.86 + temp_C))

    return esi


def tempFromEsw(esw):
    temp = (243.04 * np.log(esw / 6.1094)) / (17.625 - np.log(esw / 6.1094))

    return temp


def tempFromEsi(esi):
    temp = (273.86 * np.log(esi / 6.1121)) / (22.587 - np.log(esi / 6.1121))

    return temp


def RHFromE(e, es):
    RH = e / es

    return RH


def dewFromRH(temp_C, rh):
    # Function to calculate dewpoint given temperature and relative
    # humidity. Early soundings only measured relative humidity, not
    # dewpoint depression, and thus dewpoint must be calculated from
    # relative humidity instead of the other way around. Accuracy on this
    # calculation is very good, but note that this does not do anything to
    # change the fact that early humidity measurements from soundings could
    # often be quite inaccurate.
    #
    # General form: [dewpoint] = dewrelh(temp_C,rhum)
    #
    # Outputs:
    # dewpoint: value or vector of dewpoints (deg C)
    #
    # Inputs:
    # temp_C: value or vector of temperatures (deg C)
    # rhum: value or vector of humidity (%)

    # Coefficients from Alduchov and Eskridge 1996, accurate to relative error
    # less than 0.4 % over - 40 degC to 50 degC

    es = eswFromTemp(temp_C)
    e = (rh / 100) * es
    dewpoint_C = tempFromEsw(e)
    # a = 17.625  # dimensionless
    # b = 243.04  # deg C
    #
    # dewpoint = (b * (np.log(rh / 100) + (a * temp_C) / (b + temp_C))) / \
    #            (a - np.log(rh / 100) - (a * temp_C) / (b + temp_C))

    return dewpoint_C


def dewFromRHi(temp_C, RHi):

    esi = esiFromTemp(temp_C)
    e = (RHi / 100) * esi
    dewpoint_C = tempFromEsw(e)

    return dewpoint_C


def dewrelh(temp_C, dpd_C):
    dewpoint_C = temp_C - dpd_C  # Dewpoint is difference of temperature and dewpoint depression

    # The August-Roche-Magnus equation, accurate to within 0.4 % from -40 C to 50 C
    relative_humidity = 100 * eswFromTemp(dewpoint_C) / eswFromTemp(temp_C)

    return relative_humidity


def iceSupersatToRH(iceSupersat, temp_C):
    # Converts ice supersaturation in percent to RH in percent. Uses the
    # Improved August-Roche-Magnus saturation vapor pressure equation from:
    #  Alduchov, O.A. and R.E. Eskridge, 1996:
    #  Improved Magnus Form Approximation of Saturation Vapor Pressure.
    #  J. Appl. Meteor., 35, 601?609,
    #  https://doi.org/10.1175/1520-0450(1996)035<0601:IMFAOS>2.0.CO;2
    #  See equations 21 and 23 from above citation.
    #
    # General form: [RHw] = iceSupersatToRH(iceSupersat,temp_C)
    #
    # Output
    # RHw: relative humidity (with respect to water) in %
    #
    # Inputs:
    # iceSupersat: supersaturation with respect to ice in %
    # temp_C: temperature in deg C

    iceSupersatDecimal = iceSupersat / 100

    eswStandard = eswFromTemp(temp_C)
    esiStandard = esiFromTemp(temp_C)

    esw = esiStandard * (iceSupersatDecimal + 1)

    RHw = esw / eswStandard * 100

    return RHw


def iceSupersatToVaporExc(iceSupersatDecimal, temp_C):
    # Converts ice supersaturation in decimal to vapor density excess in
    # g/m^3. Uses the Improved August-Roche-Magnus saturation vapor pressure equation from:
    #  Alduchov, O.A. and R.E. Eskridge, 1996:
    #  Improved Magnus Form Approximation of Saturation Vapor Pressure.
    #  J. Appl. Meteor., 35, 601?609,
    #  https://doi.org/10.1175/1520-0450(1996)035<0601:IMFAOS>2.0.CO;2
    #  See equation 23 from above citation.
    #
    # General form: [vaporDensExc] = iceSupersatToVaporExc(iceSupersatDecimal,temp_C)
    #
    # Output
    # vaporDensExc: vapor density excess in g/m^3
    #
    # Inputs:
    # iceSupersatDecimal: supersaturation with respect to ice as a decimal
    # temp_C: temperature in deg C

    Rv = 461.5  # J/(kgK)
    temp_K = temp_C + 273.15

    esiStandard = esiFromTemp(temp_C)  # Saturation vapor pressure wrt ice at input temperature
    vaporPressure = iceSupersatDecimal * 100 * esiStandard
    vaporDensExc = vaporPressure / (Rv * temp_K)
    vaporDensExc = vaporDensExc * 10 ^ 3  # Convert to g/m^3 (standard unit)

    return vaporDensExc


def rhwToRhi(rhw, temp_C):
    # Converts relative humidity with respect to water to relative humidity
    # with respect to ice
    # Equation is the Improved August-Roche-Magnus approximation from:
    #  Alduchov, O.A. and R.E. Eskridge, 1996:
    #  Improved Magnus Form Approximation of Saturation Vapor Pressure.
    #  J. Appl. Meteor., 35, 601?609,
    #  https://doi.org/10.1175/1520-0450(1996)035<0601:IMFAOS>2.0.CO;2
    #  See equations 21 and 23 from above citation.
    #
    # General form: [rhi] = rhwToRhi(percent,temp_C)
    #
    # Output
    # rhi: relative humidity with respect to ice percentage
    #
    # Input
    # rhw: relative humidity with respect to water percentage
    # temp_C: temperature in Celsius

    rhw = rhw / 100

    esw = eswFromTemp(temp_C)
    esi = esiFromTemp(temp_C)

    ew = esw * rhw
    dewpt_C = tempFromEsw(ew)
    ei = esiFromTemp(dewpt_C)
    rhi = ei / esi
    rhi = rhi * 100

    return rhi


def rhiToRHw(rhi, temp_C):
    rhi = rhi / 100

    esw = eswFromTemp(temp_C)
    esi = esiFromTemp(temp_C)

    ei = esi * rhi
    dewpt_C = tempFromEsi(ei)
    ew = eswFromTemp(dewpt_C)
    rhw = ew / esw
    rhw = rhw * 100

    return rhw


def rhwToVaporExc(rhw, temp_C):
    # Calculates vapor density for a given RH percent and T. Uses the
    # Improved August-Roche-Magnus saturation vapor pressure equation from:
    #  Alduchov, O.A. and R.E. Eskridge, 1996:
    #  Improved Magnus Form Approximation of Saturation Vapor Pressure.
    #  J. Appl. Meteor., 35, 601?609,
    #  https://doi.org/10.1175/1520-0450(1996)035<0601:IMFAOS>2.0.CO;2
    #
    # General form: [rhoDiff] = rhwToVaporExc(rhw,temp_C)
    #
    # Output
    # rhoDiff: vapor density excess in g/m^3
    #
    # Input
    # rhw: Relative humidity (with respect to water) in %
    # temp_C: Temperature in Celsius

    thwDecimal = rhw / 100  # Equations use decimal instead of percent
    Rv = 461.5  # J/(kgK)
    tempK = temp_C + 273.15  # Convert temperatures to Kelvin

    eswStandard = eswFromTemp(temp_C)
    esiStandard = esiFromTemp(temp_C)
    eswStandard = eswStandard * 100
    esiStandard = esiStandard * 100
    eswPercent = eswStandard * thwDecimal

    rhow = eswPercent / (Rv * tempK)
    rhoi = esiStandard / (Rv * tempK)

    rhoDiff = rhow - rhoi
    rhoDiff = rhoDiff * 10 ^ 3

    return rhoDiff


def RHi_from_T_Td(temp_C, dewpoint_C):
    esi = esiFromTemp(temp_C)
    ei = esiFromTemp(dewpoint_C)
    rhi = ei / esi
    rhi = rhi * 100
    return rhi


def FkCalc(temp_C, iceFlag=False):
    """
    Function to calculate Fk given a certain temperatures
    :param temp_C: temperature in degC
    :param iceFlag: True if dealing with ice
    :return: Fk
    """
    # convert temperature to K
    temp_K = temp_C + 273.15

    # get correct latent heat constant
    L = constants.Lv
    if iceFlag:
        L = constants.Ls

    # calculate Fk
    Fk = (L ** 2) / (constants.thermalCond(temp_K) * constants.Rv * (temp_K ** 2))

    return Fk


def FdCalc(temp_C, es, p_kpa=100):
    """
    Function to calculate Fd given a certain temperature and vapor pressure
    :param temp_C: temperature in degC
    :param es: saturation vapor pressure
    :param p_kpa: environmental pressure in kpa
    :return: Fd
    """
    # convert temperature to K
    temp_K = temp_C + 273.15

    # calculate Fd
    Fd = (constants.Rv * temp_K) / (constants.diffCoeff(temp_K, p_kpa) * es)

    return Fd


def dropVolume(diameter_m):
    volume = (4 / 3) * np.pi * ((diameter_m / 2) ** 3)

    return volume


def dropDiameter(volume):
    diameter = 2 * ((3 * volume / (4 * np.pi)) ** (1 / 3))

    return diameter


def nearest_ind(items, value):
    diff = np.abs([item - value for item in items])
    return diff.argmin(0)
