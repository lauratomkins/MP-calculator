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
initmass_g = 1e-8
dewp_C = moist_calc.dewFromRH(temp_C, RHw) # deg C
esi = moist_calc.esiFromTemp(temp_C)  # sat. vapor pressure wrt ice
e = moist_calc.eswFromTemp(dewp_C)  # vapor pressure
shape = 'disk'
relat = 'RY1989'
initradius_m = md_calc.mass2radius(initmass_g, relat)

env = {'temp_C': temp_C,
       'e': e,
       'esi': esi,
       'p_kpa': p_kpa}

#%% Grow drop
drop_mass_g = []
drop_radius_m = []
sub_flux = []

time_step = 10 # seconds
time = np.arange(0, 1e4+time_step, time_step) # seconds

imass_g = initmass_g
iradius_m = initradius_m

for itime in time:

       # store drop information
       drop_mass_g.append(imass_g)
       drop_radius_m.append(iradius_m)

       # calculate radius
       iradius_m = md_calc.mass2radius(imass_g, relat)

       # calculate sublimation flux
       iflux_kg = calculator.sublimationFlux(env, iradius_m) # kg/s
       iflux_g = iflux_kg * 10e3
       sub_flux.append(iflux_g)

       # update drop mass
       imass_g = imass_g + (iflux_g * time_step)


#%% plot

drop_mass_ug = np.array(drop_mass_g) * 1e6

fig = plt.figure(figsize=(6,6))
ax1 = plt.gca()
ax1.plot(drop_mass_ug, time, label='RH$_{ice}$ = 116% \n RH$_{water}$ = 100%', linewidth=3)

ax1.set_xlabel('Mass [micrograms]')
ax1.set_ylabel('Time [seconds]')
ax1.set_yscale('log')
ax1.set_xscale('log')
ax1.set_ylim([1e1, 1e4])
ax1.set_xlim([1e-2, 2e2])
ax1.grid()
#ax1.set_title('Riming sensitivity', loc='right')
ax1.legend()
plt.show()
