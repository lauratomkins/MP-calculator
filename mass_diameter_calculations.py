import numpy as np

def mass2radius(mass_g, flag):

    if flag in ('RY1989', 'Rogers', 'Yau'):
        radius_cm = (mass_g / 3.8e-3)**(1/2)  # mass in g to radius in cm (Rogers and Yau - disk)
    elif flag in ('Yang2000'):
        radius_cm = 0.5 * (mass_g / (8.3e-3)) ** (1 / 2.449)  # Yang (2000) hexagonal plate
    elif flag in ('Mitchell1990'):
        mass_mg = mass_g * 1e3
        radius_mm = 0.5 * (mass_mg/(0.028))**(1/2.5) # Mitchell et al. (1990) hexagonal plate
        radius_cm = radius_mm / 10
    else:
        raise KeyError
    return radius_cm
