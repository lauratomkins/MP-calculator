import calculator
import moisture_calculations
import model_calculations
import constants
import numpy as np
import matplotlib.pyplot as plt
from metpy.units import units

#%% environment:
temp_C = -10.0 # deg C
RHw = 100.0
p_kpa = 60.0
dewp_C = moisture_calculations.dewFromRH(temp_C, RHw) # deg C
esi = moisture_calculations.esiFromTemp(temp_C)  # sat. vapor pressure wrt ice
e = moisture_calculations.eswFromTemp(dewp_C)  # vapor pressure
lwc = 1.0 # g m-3
rho_i = 0.9 # [g cm -3]
rho_w = 1.0 # [g cm -3]

drop_diameter_mm = np.array([0.5, 1, 2])
drop_diameter_cm = drop_diameter_mm / 10
drop_mass_g = (np.pi * rho_w / 6) * drop_diameter_cm**3
fallspeed_cms = np.array([233, 450, 771])
fallspeed_ms = fallspeed_cms/100
d_scwater_um = np.array([20,12])
eff = np.array([[0.96, 0.93, 0.91],[0.79, 0.79, 0.76]])

grid = {'spacing': 1000}

env = {'temp_C': temp_C,
       'e': e,
       'esi': esi,
       'p_kpa': p_kpa,
       'rho_i': rho_i,
       'lwc': lwc
       }

res = np.empty_like(eff)
for idrop, diameter in np.ndenumerate(drop_diameter_mm):
       fall_speed = fallspeed_ms[idrop]
       for icloud, c_diameter in np.ndenumerate(d_scwater_um):
              env['coll_eff'] = eff[icloud, idrop]
              res[icloud, idrop] = calculator.rimingFlux(env, diameter/2, fall_speed, grid, johnson_flag=True) # g/s

#%%

drop_mass_kg = []
drop_radius_mm = []

iradius_um = 500 # sphere
iradius_mm = iradius_um/1e3
iradius_cm = iradius_mm / 10

imass_kg = moisture_calculations.dropVolume(iradius_cm*2/100)*constants.rhoI # [kg]


time_step = 1 # seconds
time = np.arange(0, 1e4+time_step, time_step) # seconds

for itime in time:

       # store
       drop_mass_kg.append(imass_kg)
       drop_radius_mm.append(iradius_mm)

       # calculate
       iflux = calculator.rimingFlux(env, iradius_mm, fall_speed_ms_corrected=-updraft, grid=grid) # kg/s
       imass_kg = imass_kg + (iflux * time_step)

       # update
       iradius_m = moisture_calculations.dropDiameter(imass_kg/constants.rhoI)/2
       iradius_mm = iradius_m*(1e3)


#%% plot

drop_mass_ug = np.array(drop_mass_kg) * 1e9
drop_mass_mg = np.array(drop_mass_kg) * 1e6

fig = plt.figure(figsize=(6,6))
ax1 = plt.gca()
ax1.plot(drop_radius_mm, time/60, linewidth=3) #label='RH$_{ice}$ = 116% \n RH$_{water}$ = 100%', )

ax1.set_xlabel('Radius [millimeters]')
ax1.set_ylabel('Time [seconds]')
#ax1.set_yscale('log')
#ax1.set_xscale('log')
ax1.set_ylim([0, 30])
#ax1.set_xlim([1e-1, 1e2])
ax1.grid()
#ax1.set_title('Riming sensitivity', loc='right')
#ax1.legend()
plt.show()

#%% scratch
#imass_g = 1e-8
#imass_kg = imass_g / 1e3 # kg
#imass_mg = imass_g * 1e3
#iradius_cm = np.sqrt(imass_g / 3.8e-3) # mass in g to radius in cm (Rogers and Yau - disk)
#imass_g = (8.3e-3) * (2*iradius_cm)**2.449
#imass_kg = imass_g / 1e3
#iradius_cm = 0.5 * (imass_g/(8.3e-3))**(1/2.449) # Yang (2000) hexagonal plate
#iradius_mm = iradius_cm * 10 # cm to mm
#iradius_mm = 0.5 * (imass_mg/(0.028))**(1/2.5) # Mitchell et al. (1990) hexagonal plate

#iradius_cm = np.sqrt(imass_g / 3.8e-3) # Houghton (1985) disk via Rogers and Yau
#iradius_cm = 0.5 * (imass_g / (8.3e-3)) ** (1 / 2.449) # Yang (2000) hexagonal plate
# iradius_mm = 0.5 * (imass_mg / (0.028)) ** (1 / 2.5) # Mitchell et al. (1990) hexagonal plate
