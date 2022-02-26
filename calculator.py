import numpy as np
import moisture_calculations
import model_calculations
import constants
from metpy.units import units


def calc_v0(time_step, grid, starting_dsd, environment):
    # set starting values
    iheight_m = grid['top']
    idiameter_mm = starting_dsd['diameter']
    idropmass_kg = starting_dsd['mass']

    # set empty arrays for storing variables
    model_times = np.arange(0 + time_step, 10000 + time_step, time_step)
    col_mass_flux = np.zeros_like(model_times, dtype=float)
    col_drop_mass = np.zeros_like(model_times, dtype=float)
    col_drop_dia = np.zeros_like(model_times, dtype=float)
    col_latent_e = np.zeros_like(model_times, dtype=float)
    col_iwc = np.zeros_like(model_times, dtype=float)
    drop_height = np.zeros_like(model_times, dtype=float)

    # while iter_height > 0:
    # for idx, iheight in enumerate(column_heights[::-1]):
    for time_idx, itime in np.ndenumerate(model_times):
        print(str(itime) + ' seconds')
        print(str(iheight) + ' meters')

        if iheight == grid['bottom']:
            print("Reached bottom of grid; exiting loop")

            output = {
                'calc_time': model_times,
                'mass_flux': col_mass_flux,
                'drop_mass_kg': col_drop_mass,
                'drop_dia': col_drop_dia,
                'latent_e': col_latent_e,
                'iwc': col_iwc}

            return output

            break

        height_idx = moisture_calculations.nearest_ind(grid['grid_mids'],
                                                       iheight)  # np.where(grid_tops == iheight)[0][0]

        # Calculate fall speed
        fall_speed_ms = starting_dsd['fall_speed'] # m/s
        fall_speed_ms_corrected = model_calculations.fallspeedCorrection(fall_speed_ms,
                                                                      environment['air_density'][height_idx])

        # Calculate sublimation mass flux
        sublimation_flux = sublimationFlux(environment, idiameter_mm/2)

        # Calculate riming mass flux
        riming_flux = rimingFlux(environment, idiameter_mm/2, fall_speed_ms_corrected, grid)

        # update particle
        total_flux = (sublimation_flux + riming_flux) * time_step # [kg/s] * [s] = [kg]
        new_dropmass_kg = idropmass_kg + total_flux  # [kg]
        if new_dropmass_kg <= 0:
            print("Drop mass <= 0; exiting loop")

            output = {
                'calc_time': model_times,
                'mass_flux': col_mass_flux,
                'drop_mass_kg': col_drop_mass,
                'drop_dia': col_drop_dia,
                'latent_e': col_latent_e,
                'iwc': col_iwc}

            return output

            break
        #new_diameter = moisture_calculations.dropDiameter(idropmass / constants.rhoI) / 1e-3  # [mm]
        new_dropmass_g = new_dropmass_kg * 1000
        new_radius_cm = np.sqrt(new_dropmass_g / 3.8e-3) # plate crystals kg to mm
        new_diameter_mm = new_radius_cm * 2 * 10 # mm

        # store variables
        col_mass_flux[time_idx] = total_flux
        col_drop_mass[time_idx] = new_dropmass_kg
        col_drop_dia[time_idx] = new_diameter_mm
        col_latent_e[time_idx] = sublimation_flux * time_step * constants.Ls
        col_iwc[time_idx] = new_dropmass_kg / 1000e3
        drop_height[time_idx] = iheight

        # update for next iteration
        idropmass_kg = new_dropmass_kg
        idiameter_mm = new_diameter_mm
        iheight = iheight - time_step * fall_speed_ms_corrected

        if iheight <= 0:
            print("Height <= 0; exiting loop")

            output = {
                'calc_time': model_times,
                'mass_flux': col_mass_flux,
                'drop_mass_kg': col_drop_mass,
                'drop_dia': col_drop_dia,
                'latent_e': col_latent_e,
                'iwc': col_iwc}

            return output

            break

    output = {
        'calc_time': model_times,
        'mass_flux': col_mass_flux,
        'drop_mass_kg': col_drop_mass,
        'drop_dia': col_drop_dia,
        'latent_e': col_latent_e,
        'iwc': col_iwc}

    return output

def sublimationFlux(env, drop_radius_m, type='disk'):
    """

    :param env:
    :param drop_radius_m:
    :param type:
    :return:
    """
    # input array, output array
    # input value, output value

    supersat = (env['e'] / env['esi']) - 1
    Fk = moisture_calculations.FkCalc(env['temp_C'], iceFlag=True)  # [m s kg-1]
    Fd = moisture_calculations.FdCalc(env['temp_C'], 100*env['esi'], p_kpa=env['p_kpa'])  # [m s kg-1]

    # Calculate sublimation mass flux
    if type in ('disk'):
        capacitance = 2 * (drop_radius_m) / np.pi  # disk 2r/pi sphere c=r [m]
    elif type in ('sphere'):
        capacitance = drop_radius_m

    # Houze (2003) textbook Ch. 3 and Rogers and Yau (1989)
    sublimation_flux = (4 * np.pi * capacitance * supersat) / (Fk + Fd) # [kg/s]

    return sublimation_flux

def rimingFlux(env, drop_radius_m, fall_speed_ms_corrected, grid=None, johnson_flag=False):

    drop_radius_mm = drop_radius_m * 1e3

    if johnson_flag:
        rho_w = 1.0 # [g cm-3]
        r = 1 # axis ratio
        drop_diameter_cm = (drop_radius_mm * 2)/10
        drop_mass_g = (np.pi * rho_w / 6) * drop_diameter_cm**3
        A_cm = (np.pi/4) * (6 * drop_mass_g / (np.pi * r * env['rho_i'])) ** (2/3) # [cm2]
        A_m = A_cm / (100**2) # [m2]
    else:
        A_m = np.pi * ((drop_radius_mm + (env['scwater_d_mm'] / 2)) * 1e-3) ** 2 # cross sectional area [m2]

    if 'lwc' in env.keys():
        lwc = env['lwc'] # [g m3]
    else:
        lwc = (env['scwater_n'] * moisture_calculations.dropVolume(env['scwater_d_mm'] * 1e-3) * constants.rhoW) \
              / (grid['spacing'] * 1000 ** 2) # ndrops * mass of 1 drop / grid volume [g / m3]

    # Johnson (1987) and Houze (2003) textbook Ch. 3
    riming_flux = A_m * abs(fall_speed_ms_corrected - 0) * env['coll_eff'] * lwc # [g/s]


    return riming_flux

def createGrid(top=10000, bottom=0, spacing=1000):

    # All values in meters

    grid_tops, grid_bots, grid_mids = model_calculations.model_grid(bottom, top, spacing) # m

    grid = {
        'top': top,
        'bottom': bottom,
        'spacing': spacing,
        'grid_tops': grid_tops,
        'grid_mids': grid_mids,
        'grid_bots': grid_bots}

    return grid

def createDSD(ndrops, diameter_mm, fall_speed):

    #drop_mass_kg = ndrops * moisture_calculations.dropVolume(diameter * 1e-3) * constants.rhoI  # [kg]
    diameter_cm = diameter_mm / 10
    radius_cm = diameter_cm / 2
    drop_mass_g = (3.8e-3 * (radius_cm ** 2)) # Houghton (1985) via Rogers and Yau # cm to g
    drop_mass_kg = drop_mass_g / 1000

    starting_dsd = {
        'diameter': diameter_mm,
        'mass': drop_mass_kg,
        'fall_speed': fall_speed}

    return starting_dsd

def createEnv(grid, tempLims=[0,-30], constantRHi=100, sc_layers=False, isclayer=7, nsclayer=100, dsclayer=0.01):
    grid_mids = grid['grid_mids']

    # Create temperature, dewpoint, RHi profile
    temp_C = np.linspace(tempLims[0], tempLims[1], len(grid_mids))
    RHi = np.full(np.shape(temp_C), constantRHi)
    dewp_C = moisture_calculations.dewFromRHi(temp_C, RHi)
    esi = moisture_calculations.esiFromTemp(temp_C)  # sat. vapor pressure wrt ice
    e = moisture_calculations.eswFromTemp(dewp_C)  # vapor pressure
    es = moisture_calculations.eswFromTemp(temp_C)  # sat. vapor pressure wrt liquid
    RHw = e / es * 100
    rho_mids = model_calculations.ZtoAirDensity(grid_mids)  # [kg/m3]

    # Defining supercooled liquid water layer for riming
    n_scwater = np.zeros_like(RHi, dtype=float)
    d_scwater_mm = np.zeros_like(RHi, dtype=float)
    collection_efficiency = 0
    if sc_layers:
        n_cm = nsclayer  # droplets/cm3
        n_m = n_cm * (100 ** 3) # droplets/m3
        n_scwater[isclayer] = n_m * (grid['spacing'] ** 3)
        d_scwater_mm[isclayer] = dsclayer  # [mm] (10 micrometers)
        collection_efficiency = 1

    environment = {
        'temp_C': temp_C,
        'dtemp_C': dewp_C,
        'RHw': RHw,
        'RHi': RHi,
        'e': e,
        'es': es,
        'esi': esi,
        'scwater_n': n_scwater,
        'scwater_d_mm': d_scwater_mm,
        'coll_eff': collection_efficiency,
        'air_density': rho_mids,
        'height': grid_mids}

    return environment