import calculator
import moisture_calculations as moist_calc
import mass_diameter_calculations as md_calc
import model_calculations
import constants
import numpy as np
import matplotlib.pyplot as plt

#%% vapor deposition environment:
# temperature in degC
temp_C = -15
# RH wrt ice (testing 3 values)
RHi = np.array([110, 120, 130])
# pressure in kPa
p_kpa = 80
# initial mass of crystal in g
initmass_g = 1e-8
# calculate dewpoint from temperature and RHice in degC
dewp_C = moist_calc.dewFromRHi(temp_C, RHi)
# sat. vapor pressure wrt ice from temperature
esi = moist_calc.esiFromTemp(temp_C)
# vapor pressure from dewpoint
e = moist_calc.eswFromTemp(dewp_C)
# mass-radius relationship to use
relat = 'RY1989'
# shape of particle (optional - only some relationships need this)
shape = 'disk'
# initial radius calculated from mass
initradius_m = md_calc.mass2radius(initmass_g, relat)
# plotting colors for deposition lines
dep_colors = ['#ffc1b9', '#ff948c', '#ff6361']

# initialize empty dictionary to store information
dep_env = dict()

# loop through vapor deposition environments and save
for ind, vpres in enumerate(e):
   env_str = 'env' + str(RHi[ind])
   dep_env[env_str] = {'temp_C': temp_C,
                       'e': vpres,
                       'esi': esi,
                       'p_kpa': p_kpa,
                       'RHi': str(RHi[ind]),
                       'color': dep_colors[ind]}

#%% riming environment:
# temperature in degC
temp_C = -15
# constant RHice (no vapor deposition)
RHi = 100
# environment pressure in kPa
p_kpa = 80
# inital mass in g
initmass_g = 1e-8
# dewpoint
dewp_C = moist_calc.dewFromRHi(temp_C, RHi)
# sat. vapor pressure wrt ice
esi = moist_calc.esiFromTemp(temp_C)
# vapor pressure
e = moist_calc.eswFromTemp(dewp_C)
# mass-radius relationship to use
relat = 'RY1989'
# shape of particle (optional - only some relationships need this)
shape = 'disk'
# inital radius in meters
initradius_m = md_calc.mass2radius(initmass_g, relat)
# liquid water equivalent (3 values)
lwc = np.array([0.5, 1.0, 2.0]) # g m-3
# density of ice
rho_i = 0.9 # [g cm -3]
# plotting colors
rim_colors = ['#00c6e9','#008ec3','#075891']

# initialize empty dictionary
rim_env = dict()

# loop through riming environmetns and save in dictionary
for lw_ind, content in enumerate(lwc):
   env_str = 'env' + str(content)
   rim_env[env_str] = {'temp_C': temp_C,
                       'e': e,
                       'esi': esi,
                       'p_kpa': p_kpa,
                       'lwc': content,
                       'rho_i': rho_i,
                       'coll_eff': 0.95,
                       'RHi': str(RHi),
                       'color': rim_colors[lw_ind]}

#%% Grow drop by vapor deposition
dep_drop = dict()
# Loop through environments
for rh_env in dep_env.keys():

    # isolate environment
   env = dep_env[rh_env]

   # set empty arrays
   drop_mass_g = []
   drop_radius_m = []
   sub_flux = []

   # set time steps
   time_step = 1 # seconds
   time = np.arange(0, 1e4+time_step, time_step) # seconds

   # initialize drop
   imass_g = initmass_g
   iradius_m = initradius_m

   # run through time
   for itime in time:

          # store drop information
          drop_mass_g.append(imass_g)
          drop_radius_m.append(iradius_m)

          # calculate radius
          iradius_m = md_calc.mass2radius(imass_g, relat)

          # calculate sublimation flux
          iflux_kg = calculator.sublimationFlux(env, iradius_m, type=shape) # kg/s
          iflux_g = iflux_kg * 1e3
          sub_flux.append(iflux_g)

          # update drop mass
          imass_g = imass_g + (iflux_g * time_step)

   # save information
   dep_drop[rh_env] = {'drop_mass_g': np.array(drop_mass_g),
                        'drop_radius_m': np.array(drop_radius_m),
                        'sub_flux': np.array(sub_flux),
                        'time': np.array(time)}

#%% Grow drop by riming
fall_speed_ms = 1 # m/s
rim_drop = dict()
# Loop through environments
for lwc_env in rim_env.keys():

   env = rim_env[lwc_env]

   # set empty arrays
   drop_mass_g = []
   drop_radius_m = []
   sub_flux = []

   # set time steps
   time_step = 1 # seconds
   time = np.arange(0, 1e4+time_step, time_step) # seconds

   # initialize drop
   imass_g = initmass_g
   iradius_m = initradius_m

   # run through time
   for itime in time:

          # store drop information
          drop_mass_g.append(imass_g)
          drop_radius_m.append(iradius_m)

          # calculate radius
          iradius_m = md_calc.mass2radius(imass_g, relat)

          # calculate sublimation flux
          iflux_g = calculator.rimingFlux(env, iradius_m, fall_speed_ms,
                                           grid=None, johnson_flag=True) # kg/s
          sub_flux.append(iflux_g)

          # update drop mass
          imass_g = imass_g + (iflux_g * time_step)

   # save information
   rim_drop[lwc_env] = {'drop_mass_g': np.array(drop_mass_g),
                        'drop_radius_m': np.array(drop_radius_m),
                        'sub_flux': np.array(sub_flux),
                        'time': np.array(time)}

#%% plot

fig = plt.figure(figsize=(6,6))
ax1 = plt.gca()
for key in dep_drop.keys():
   ax1.plot(dep_drop[key]['time'], dep_drop[key]['drop_mass_g'] * 1e6,
            label='RH$_{{ice}}$ = {0}%'.format(dep_env[key]['RHi']),
            color=dep_env[key]['color'], linewidth=3)

for key in rim_drop.keys():
   ax1.plot(rim_drop[key]['time'], rim_drop[key]['drop_mass_g'] * 1e6,
            label='LWC = {0} g m$^{{-3}}$'.format(str(rim_env[key]['lwc'])),
            color=rim_env[key]['color'], linewidth=3)

ax1.set_ylabel('Mass [micrograms]')
ax1.set_xlabel('Time [seconds]')
ax1.set_yscale('log')
ax1.set_xscale('log')
ax1.set_xlim([1e0, 1e4])
ax1.set_ylim([1e-2, 1e7])
ax1.grid(linestyle='--')
#ax1.set_title('Riming sensitivity', loc='right')
ax1.legend()
#plt.savefig('Q:\\My Drive\\phd\\exams\\oral_exam\\figs\\mp_calc.png', dpi=600, bbox_inches='tight')
plt.show()
