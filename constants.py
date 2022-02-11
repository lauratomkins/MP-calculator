import numpy as np

# Latent Heat constants (ametsoc glossary)
Lv = 2.501e6  # J/kg
Lf = 3.337e5  # J/kg
Ls = 2.834e6  # J/kg

# Gas constant for water vapor (ametsoc glossary)
Rv = 461.5  # J/(kgK)

# density of water (liquid and ice)
rhoW = 999.8  # kg/m3
rhoI = 917  # kg/m3


# Thermal conductivity of air (Ka [J m-1 s-1 K-1]) Rogers and Yau p. 103
def thermalCond(temp_K):
    temps = np.array([-40, -30, -20, -10, 0, 10, 20, 30]) # degC
    tempsK = temps + 273.15
    Ka = np.array([2.07e-2, 2.16e-2, 2.24e-2, 2.32e-2, 2.4e-2, 2.48e-2, 2.55e-2, 2.63e-2])
    mk, bk = np.polyfit(tempsK, Ka, 1)

    Ka = (mk * temp_K) + bk
    return Ka


# Diffusion coefficient (Dv [m2 s-1]) Rogers and Yau p. 103
def diffCoeff(temp_K, p_kpa=100):
    temps = np.array([-40, -30, -20, -10, 0, 10, 20, 30]) # degC
    tempsK = temps + 273.15
    Dv = np.array([1.62e-5, 1.76e-5, 1.91e-5, 2.06e-5, 2.21e-5, 2.36e-5, 2.52e-5, 2.69e-5]) * (100 / p_kpa)
    md, bd = np.polyfit(tempsK, Dv, 1)

    Dv = (md * temp_K) + bd
    return Dv
