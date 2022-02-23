import calculator
import moisture_calculations as moist_calc
import mass_diameter_calculations as md_calc
import model_calculations
import constants
import numpy as np
import matplotlib.pyplot as plt

#%% environment:
temp_C = -15 # deg C
RHw = 100
p_kpa = 80
imass_g = 1e-8
dewp_C = moist_calc.dewFromRH(temp_C, RHw) # deg C
esi = moist_calc.esiFromTemp(temp_C)  # sat. vapor pressure wrt ice
e = moist_calc.eswFromTemp(dewp_C)  # vapor pressure
shape = 'disk'
relat = 'Mitchell1990'

env = {'temp_C': temp_C,
       'e': e,
       'esi': esi,
       'p_kpa': p_kpa}

drop_mass_kg = []
drop_radius_mm = []
sub_flux = []

iradius_cm = md_calc.mass2radius(imass_g, relat)
iradius_mm = iradius_cm * 10
imass_kg = imass_g / 1e3 # kg

time_step = 100 # seconds
time = np.arange(0, 1e4+time_step, time_step) # seconds

for itime in time:

       drop_mass_kg.append(imass_kg)
       drop_radius_mm.append(iradius_mm)

       iradius_cm = md_calc.mass2radius(imass_g, relat)
       iradius_mm = iradius_cm * 10

       iflux = calculator.sublimationFlux(env, iradius_mm) # kg/s

       sub_flux.append(iflux)
       imass_kg = imass_kg + (iflux * time_step)
       imass_g = imass_kg * 1e3


#%% plot

drop_mass_ug = np.array(drop_mass_kg) * 1e9

fig = plt.figure(figsize=(6,6))
ax1 = plt.gca()
ax1.plot(drop_mass_ug, time, label='RH$_{ice}$ = 116% \n RH$_{water}$ = 100%', linewidth=3)

ax1.set_xlabel('Mass [micrograms]')
ax1.set_ylabel('Time [seconds]')
ax1.set_yscale('log')
ax1.set_xscale('log')
ax1.set_ylim([1e0, 1e4])
ax1.set_xlim([1e-2, 1e2])
ax1.grid()
#ax1.set_title('Riming sensitivity', loc='right')
ax1.legend()
plt.show()
