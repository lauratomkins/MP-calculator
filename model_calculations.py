import numpy as np
#import moisture_calculations

def model_grid(bottom, top, spacing):

    grid_tops = np.arange(bottom + spacing, top + spacing, spacing)
    grid_bots = np.arange(bottom, top, spacing)
    grid_mids = np.arange(bottom + 0.5 * spacing, top + 0.5 * spacing, spacing)
    #rho_mids = moisture_calculations.ZtoAirDensity(grid_mids) # [kg/m3]

    return grid_tops, grid_bots, grid_mids

def fall_speed_from_d(diameter): # FOR GRAUPEL/RAIN NOT SNOW!!

    fall_speed = 1.6 * (diameter ** 0.46)

    return fall_speed

def fallspeedCorrection(V, rho):
    # From D. Kingsmill code
    V_corrected = V * (1.275 / rho) ** 0.4

    return V_corrected


def ZtoAirDensity(Z):
    # From D. Kingsmill code. Z in meters, rho in kg/m3
    rho = 1.225 * np.exp(Z * -1e-4)

    return rho