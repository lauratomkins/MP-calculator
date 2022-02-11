import calculator
import moisture_calculations
import model_calculations
import constants
import numpy as np
import matplotlib.pyplot as plt

#%% environment:
temp_C = -15 # deg C
RHw = 100
p_kpa = 80
dewp_C = moisture_calculations.dewFromRH(temp_C, RHw) # deg C
esi = moisture_calculations.esiFromTemp(temp_C)  # sat. vapor pressure wrt ice
e = moisture_calculations.eswFromTemp(dewp_C)  # vapor pressure

env = {'temp_C': temp_C,
       'e': e,
       'esi': esi,
       'p_kpa': p_kpa}

drop_mass_kg = []
drop_radius_mm = []

iradius_mm = 1
iradius_cm = iradius_mm / 10
#imass_g = 1e-8
#imass_kg = imass_g / 1e3 # kg
#imass_mg = imass_g * 1e3
#iradius_cm = np.sqrt(imass_g / 3.8e-3) # mass in g to radius in cm (Rogers and Yau - disk)
imass_g = (8.3e-3) * (2*iradius_cm)**2.449
imass_kg = imass_g / 1e3
#iradius_cm = 0.5 * (imass_g/(8.3e-3))**(1/2.449) # Yang (2000) hexagonal plate
#iradius_mm = iradius_cm * 10 # cm to mm
#iradius_mm = 0.5 * (imass_mg/(0.028))**(1/2.5) # Mitchell et al. (1990) hexagonal plate

time_step = 1 # seconds
time = np.arange(0, 1e4+time_step, time_step) # seconds

for itime in time:

       drop_mass_kg.append(imass_kg)
       drop_radius_mm.append(iradius_mm)

       #iradius_cm = np.sqrt(imass_g / 3.8e-3) # Houghton (1985) disk via Rogers and Yau
       iradius_cm = 0.5 * (imass_g / (8.3e-3)) ** (1 / 2.449) # Yang (2000) hexagonal plate
       iradius_mm = iradius_cm * 10
       #iradius_mm = 0.5 * (imass_mg / (0.028)) ** (1 / 2.5) # Mitchell et al. (1990) hexagonal plate
       iflux = calculator.sublimationFlux(env, iradius_mm) # kg/s
       imass_kg = imass_kg + (iflux * time_step)
       imass_g = imass_kg * 1e3
       imass_mg = imass_g * 1e3


#%% plot

drop_mass_ug = np.array(drop_mass_kg) * 1e9
drop_mass_mg = np.array(drop_mass_kg) * 1e6

fig = plt.figure(figsize=(6,6))
ax1 = plt.gca()
ax1.plot(drop_mass_mg, time, label='RH$_{ice}$ = 116% \n RH$_{water}$ = 100%', linewidth=3)

ax1.set_xlabel('Mass [milligrams]')
ax1.set_ylabel('Time [seconds]')
ax1.set_yscale('log')
ax1.set_xscale('log')
ax1.set_ylim([1e0, 1e5])
ax1.set_xlim([1e-1, 1e2])
ax1.grid()
#ax1.set_title('Riming sensitivity', loc='right')
ax1.legend()
plt.show()