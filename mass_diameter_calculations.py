"""
Function(s) for mass-diameter relationships
"""
import numpy as np

def mass2radius(mass_g, relat_flag='RY1989'):
    """
    Function to convert droplet mass to droplet radius
    :param mass_g: mass of drops in grams
    :param relat_flag: string of relationshop to use
    :return: radius in meters
    """

    if relat_flag in ('RY1989', 'Rogers', 'Yau'):
        radius_cm = (mass_g / 3.8e-3)**(1/2)  # mass in g to radius in cm (Rogers and Yau - disk)
    elif relat_flag in ('Yang2000'):
        radius_cm = 0.5 * (mass_g / (8.3e-3)) ** (1 / 2.449)  # Yang (2000) hexagonal plate
    elif relat_flag in ('Mitchell1990'):
        mass_mg = mass_g * 1e3
        radius_mm = 0.5 * (mass_mg/(0.028))**(1/2.5) # Mitchell et al. (1990) hexagonal plate
        radius_cm = radius_mm / 10
    else:
        raise KeyError

    radius_m = radius_cm * 1e-2

    return radius_m
