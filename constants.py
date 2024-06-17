"""
constants used in calculations
"""
import numpy as np

# Latent Heat constants (ametsoc glossary)
Lv = 2.501e6  # J/kg
Lf = 3.337e5  # J/kg
Ls = 2.834e6  # J/kg

# Gas constant for water vapor (ametsoc glossary)
Rv = 461.5  # J/(kgK)
Rd = 287    # J/(kgK)

# density of water (liquid and ice)
rhoW = 999.8  # kg/m3
rhoI = 917  # kg/m3


# Thermal conductivity of air (Ka [J m-1 s-1 K-1]) Rogers and Yau p. 103
def thermalCond(temp_K_in):
    """
    Returns thermal conductivity of air (Ka) for a given temperature
    :param temp_K_in: Temperature in Kelvin
    :return: Ka_out: Ka for specified temperature
    """
    # temperatures and Ka for model
    temps = np.array([-40, -30, -20, -10, 0, 10, 20, 30]) # degC
    tempsK = temps + 273.15 # Kelvin
    Ka = np.array([2.07e-2, 2.16e-2, 2.24e-2, 2.32e-2, 2.4e-2, 2.48e-2, 2.55e-2, 2.63e-2])

    # fit line
    mk, bk = np.polyfit(tempsK, Ka, 1)

    # calculate Ka
    Ka_out = (mk * temp_K_in) + bk

    return Ka_out


# Diffusion coefficient (Dv [m2 s-1]) Rogers and Yau p. 103
def diffCoeff(temp_K_in, p_kpa=100):
    """
    Returns diffusion coefficient (Dv) for a given temperature
    :param temp_K: temperature in Kelvin
    :param p_kpa: pressure in kpa
    :return: Dv
    """
    # temperatures and Dv for model
    temps = np.array([-40, -30, -20, -10, 0, 10, 20, 30]) # degC
    tempsK = temps + 273.15 # K
    Dv = np.array([1.62e-5, 1.76e-5, 1.91e-5, 2.06e-5, 2.21e-5, 2.36e-5, 2.52e-5, 2.69e-5]) * (100 / p_kpa)

    # fit line
    md, bd = np.polyfit(tempsK, Dv, 1)

    # calculate Dv
    Dv_out = (md * temp_K_in) + bd

    return Dv_out
