import numpy as np
import moisture_calculations
import model_calculations
import constants


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

        # Find supersaturation, Fk, and Fd at height
        supersat = (environment['e'][height_idx] / environment['esi'][height_idx]) - 1
        Fk = moisture_calculations.FkCalc(environment['temp_C'][height_idx], iceFlag=True)  # [m s kg-1]
        Fd = moisture_calculations.FdCalc(environment['temp_C'][height_idx], environment['esi'][height_idx])  # [m s kg-1]

        # Calculate fall speed
        fall_speed_ms = starting_dsd['fall_speed'] # m/s
        fall_speed_ms_corrected = model_calculations.fallspeedCorrection(fall_speed_ms,
                                                                      environment['air_density'][height_idx])

        # Calculate sublimation mass flux
        capacitance = 2 * ((idiameter_mm / 2) * 1e-3) / np.pi # disk 2r/pi sphere c=r
        sublimation_flux = (4 * np.pi * capacitance * supersat) / (Fk + Fd)  # dm/dt [kg s-1]

        # Calculate riming mass flux
        A_m = np.pi * (((idiameter_mm / 2) + (environment['scwater']['d_mm'][height_idx] / 2)) * 1e-3) ** 2
        riming_flux = A_m * abs(fall_speed_ms_corrected - 0) * environment['scwater']['coll_eff'] * \
                      (environment['scwater']['n'][height_idx] * moisture_calculations.dropVolume(
                          environment['scwater']['d_mm'][height_idx] * 1e-3) * constants.rhoI) / (
                              grid['spacing'] * 1000 ** 2)

        #print(str(riming_flux))

        # update particle
        total_flux = (sublimation_flux + riming_flux) * time_step
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

def sublimationFlux(env, drop_radius_mm):

    supersat = (env['e'] / env['esi']) - 1
    Fk = moisture_calculations.FkCalc(env['temp_C'], iceFlag=True)  # [m s kg-1]
    Fd = moisture_calculations.FdCalc(env['temp_C'], env['esi'], env['p_kpa'])  # [m s kg-1]

    # Calculate sublimation mass flux
    capacitance = 2 * (drop_radius_mm) / np.pi  # disk 2r/pi sphere c=r [m]
    sublimation_flux = (4 * np.pi * capacitance * supersat) / (Fk + Fd) # [kg/s]

    return sublimation_flux

def rimingFlux(env, drop_radius_mm, fall_speed_ms_corrected):

    A_m = np.pi * ((drop_radius_mm + (env['scwater']['d_mm'] / 2)) * 1e-3) ** 2
    riming_flux = A_m * abs(fall_speed_ms_corrected - 0) * env['scwater']['coll_eff'] * \
                  (env['scwater']['n'] * moisture_calculations.dropVolume(
                      env['scwater']['d_mm'] * 1e-3) * constants.rhoI) / (
                          grid['spacing'] * 1000 ** 2)


def createGrid(top=10000, bottom=0, spacing=1000):

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
        'scwater': {
            'n': n_scwater,
            'd_mm': d_scwater_mm,
            'coll_eff': collection_efficiency},
        'air_density': rho_mids}

    return environment