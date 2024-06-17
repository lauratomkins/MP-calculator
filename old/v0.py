"""
Early versions of the calculator
"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import FixedLocator

import calculator
import moisture_calculations
import model_calculations
import constants

# %% Define grid
import plotting

grid = calculator.createGrid(10000, 0, 1000)
grid_rogers_yau = calculator.createGrid(10000, 0, 100)
time_step = 1000 # s

# %% dsds
dsd_1mm_norm = calculator.createDSD(ndrops=1, diameter=1, fall_speed=1) # 1 x 1 mm drop
dsd_1mm_slow = calculator.createDSD(ndrops=1, diameter=1, fall_speed=0.75) # 1 x 1 mm drop
dsd_1mm_fast = calculator.createDSD(ndrops=1, diameter=1, fall_speed=1.25) # 1 x 1 mm drop

dsd_2mm_norm = calculator.createDSD(ndrops=1, diameter=2, fall_speed=1) # 1 x 1 mm drop
dsd_2mm_slow = calculator.createDSD(ndrops=1, diameter=2, fall_speed=0.75) # 1 x 1 mm drop
dsd_2mm_fast = calculator.createDSD(ndrops=1, diameter=2, fall_speed=1.25) # 1 x 1 mm drop

dsd_4mm_norm = calculator.createDSD(ndrops=1, diameter=4, fall_speed=1) # 1 x 1 mm drop
dsd_4mm_slow = calculator.createDSD(ndrops=1, diameter=4, fall_speed=0.75) # 1 x 1 mm drop
dsd_4mm_fast = calculator.createDSD(ndrops=1, diameter=4, fall_speed=1.25) # 1 x 1 mm drop

dsd_rogers_yau = calculator.createDSD(ndrops=1, diameter=0.0324, fall_speed=1)

# %% environments
env_rhi101 = calculator.createEnv(grid, tempLims=[0,-30], constantRHi=101, sc_layers=False)
env_rhi102 = calculator.createEnv(grid, tempLims=[0,-30], constantRHi=102.5, sc_layers=False)
env_rhi105 = calculator.createEnv(grid, tempLims=[0,-30], constantRHi=105, sc_layers=False)
env_rhi110 = calculator.createEnv(grid, tempLims=[0,-30], constantRHi=110, sc_layers=False)
env_rhi120 = calculator.createEnv(grid, tempLims=[0,-30], constantRHi=120, sc_layers=False)
env_rhi130 = calculator.createEnv(grid, tempLims=[0,-30], constantRHi=130, sc_layers=False)

env_n100 = calculator.createEnv(grid, tempLims=[0,-30], constantRHi=100,
                                sc_layers=True, isclayer=7, nsclayer=100, dsclayer=0.01)
env_n500 = calculator.createEnv(grid, tempLims=[0,-30], constantRHi=100,
                                sc_layers=True, isclayer=7, nsclayer=500, dsclayer=0.01)
env_n1000 = calculator.createEnv(grid, tempLims=[0,-30], constantRHi=100,
                                sc_layers=True, isclayer=7, nsclayer=1000, dsclayer=0.01)

env_rogers_yau = calculator.createEnv(grid_rogers_yau, tempLims=[0, -30], constantRHi=113, sc_layers=False)
#%% riming calculator
output_n100 = calculator.calc_v0(time_step, grid, dsd_1mm_norm, env_n100)
output_n500 = calculator.calc_v0(time_step, grid, dsd_1mm_norm, env_n500)
output_n1000 = calculator.calc_v0(time_step, grid, dsd_1mm_norm, env_n1000)

# sensitivities
output_n500_r1 = output_n500
output_n500_r2 = calculator.calc_v0(time_step, grid, dsd_2mm_norm, env_n500)
output_n500_r4 = calculator.calc_v0(time_step, grid, dsd_4mm_norm, env_n500)

output_n500_norm = output_n500
output_n500_slow = calculator.calc_v0(time_step, grid, dsd_1mm_slow, env_n500)
output_n500_fast = calculator.calc_v0(time_step, grid, dsd_1mm_fast, env_n500)

output_rogers_yau = calculator.calc_v0(100, grid_rogers_yau, dsd_rogers_yau, env_rogers_yau)

# %% deposition calculator
output_rh101 = calculator.calc_v0(time_step, grid, dsd_1mm_norm, env_rhi101)
output_rh102 = calculator.calc_v0(time_step, grid, dsd_1mm_norm, env_rhi102)
output_rh105 = calculator.calc_v0(time_step, grid, dsd_1mm_norm, env_rhi105)
output_rh110 = calculator.calc_v0(time_step, grid, dsd_1mm_norm, env_rhi110)
output_rh120 = calculator.calc_v0(time_step, grid, dsd_1mm_norm, env_rhi120)
output_rh130 = calculator.calc_v0(time_step, grid, dsd_1mm_norm, env_rhi130)

output_rh110_r1 = output_rh110
output_rh110_r2 = calculator.calc_v0(time_step, grid, dsd_2mm_norm, env_rhi110)
output_rh110_r4 = calculator.calc_v0(time_step, grid, dsd_4mm_norm, env_rhi110)

output_rh110_norm = output_rh110
output_rh110_slow = calculator.calc_v0(time_step, grid, dsd_1mm_slow, env_rhi110)
output_rh110_fast = calculator.calc_v0(time_step, grid, dsd_1mm_fast, env_rhi110)

#%% get ready to plot
model_times = output_rh101['calc_time']
model_times = np.insert(output_n100['calc_time'], 0, 0)

plot_type = 'radius' #'iwc', 'radius'
xunits = 'mins' # 'seconds'

# y labels
if plot_type == 'mass':
    dat_key = 'drop_mass_kg'
    ylab = 'Mass [kg]'
elif plot_type == 'iwc':
    dat_key = 'iwc'
    ylab = 'IWC [g m$^{-3}$]'
elif plot_type == 'radius':
    dat_key = 'drop_dia'
    ylab = 'Radius [mm]'

# x labels
if xunits == 'mins':
    xlab = 'Time [minutes]'
    model_times = model_times/60
elif xunits == 'seconds':
    xlab = 'Time [seconds]'

# sensitivity labels
radii = ['0.5 mm', '1 mm', '2 mm']
speed = ['1 ms$^{-1}$', '2 ms$^{-1}$', '0.75 ms$^{-1}$']
#%% plotting dictionaries
# riming sensitivites
riming_s_sensitivity = {'norm': output_n500_norm, 'fast': output_n500_fast, 'slow': output_n500_slow}
riming_s_plot = {}
for key in riming_s_sensitivity:
    riming_s_plot[key] = {'data': riming_s_sensitivity[key][dat_key],
                          'label': speed[list(riming_s_sensitivity.keys()).index(key)]}

riming_r_sensitivity = {'r1': output_n500_r1, 'r2': output_n500_r2, 'r4': output_n500_r4}
riming_r_plot = {}
for key in riming_r_sensitivity:
    riming_r_plot[key] = {'data': riming_r_sensitivity[key][dat_key]/2,
                          'label': radii[list(riming_r_sensitivity.keys()).index(key)]}

# deposition sensitivies
deposition_s_sensitivity = {'norm': output_rh110_norm, 'fast': output_rh110_fast, 'slow': output_rh110_slow}
deposition_s_plot = {}
for key in deposition_s_sensitivity:
    deposition_s_plot[key] = {'data': deposition_s_sensitivity[key][dat_key],
                              'label': speed[list(deposition_s_sensitivity.keys()).index(key)]}

deposition_r_sensitivity = {'r1': output_rh110_r1, 'r2': output_rh110_r2, 'r4': output_rh110_r4}
deposition_r_plot = {}
for key in deposition_r_sensitivity:
    deposition_r_plot[key] = {'data': deposition_r_sensitivity[key][dat_key]/2,
                              'label': radii[list(deposition_r_sensitivity.keys()).index(key)]}

#%% Riming Radii sensitivity

fig = plt.figure(figsize=(6,6))
ax = plt.gca()
plotting.plot_curves(ax, model_times, riming_r_plot, ylab=ylab, xlab=xlab, title='Riming radii sensitivity')
#plt.show()
plt.savefig("Q:\\My Drive\\phd\\winter_storms\\seeder_feeder\\symposium\\riming_r.png", dpi=300, bbox_inches='tight')


#%% Riming Speed sensitivity

fig = plt.figure(figsize=(6,6))
ax = plt.gca()
plotting.plot_curves(ax, model_times, riming_s_plot, use_first= False, first_val=1, ylab=ylab, xlab=xlab, title='Riming speed sensitivity')
#plt.show()
plt.savefig("Q:\\My Drive\\phd\\winter_storms\\seeder_feeder\\symposium\\riming_s.png", dpi=300, bbox_inches='tight')


#%% Deposition Radii sensitivity

fig = plt.figure(figsize=(6,6))
ax = plt.gca()
plotting.plot_curves(ax, model_times, deposition_r_plot, ylab=ylab, xlab=xlab, title='Deposition radii sensitivity')
#plt.show()
plt.savefig("Q:\\My Drive\\phd\\winter_storms\\seeder_feeder\\symposium\\deposition_r.png", dpi=300, bbox_inches='tight')


#%% Deposition Speed sensitivity

fig = plt.figure(figsize=(6,6))
ax = plt.gca()
plotting.plot_curves(ax, model_times, deposition_s_plot, use_first= False, first_val=1, ylab=ylab, xlab=xlab, title='Deposition speed sensitivity')
#plt.show()
plt.savefig("Q:\\My Drive\\phd\\winter_storms\\seeder_feeder\\symposium\\deposition_s.png", dpi=300, bbox_inches='tight')


#%% Rogers and Yau sensitivity
model_times = output_rogers_yau['calc_time']
model_times = np.insert(model_times, 0, 0)

plot_dia = (output_rogers_yau['drop_dia']/2)*1000
plot_dia = np.insert(plot_dia, 0, plot_dia[0])
plot_dia = np.ma.masked_values(plot_dia, 0)

plot_mass = output_rogers_yau['drop_mass_kg']*1e9
plot_mass = np.insert(plot_mass, 0, 1e-2)
plot_mass = np.ma.masked_values(plot_mass, 0)

fig = plt.figure(figsize=(6,6))
ax1 = plt.gca()
ax1.plot(plot_mass, model_times, label='RH$_{ice}$ = 113%', linewidth=3)

ax1.set_xlabel('Mass [micrograms]')
ax1.set_ylabel('Time [seconds]')
#ax1.set_ylim([0,60])
ax1.set_yscale('log')
ax1.set_xscale('log')
ax1.set_ylim([1e1, 1e4])
ax1.set_xlim([1e-2, 1e2])
#ax1.set_xscale('log')
ax1.grid()
#ax1.set_title('Riming sensitivity', loc='right')
ax1.legend()
plt.show()

#%% calculations

model_times = np.insert(output_n100['calc_time'], 0, 0)

#d1_mass = output_d1['drop_mass_kg']; d1_mass = np.insert(d1_mass, 0, starting_drop_mass); d1_mass = np.ma.masked_values(d1_mass, 0)
#d2_mass = output_d2['drop_mass_kg']; d2_mass = np.insert(d2_mass, 0, starting_drop_mass); d2_mass = np.ma.masked_values(d2_mass, 0)
d4_mass = output_d4['drop_mass_kg']; d4_mass = np.insert(d4_mass, 0, starting_drop_mass); d4_mass = np.ma.masked_values(d4_mass, 0)

interp_times = np.arange(0, 10200, 200)

d1_interp = np.interp(interp_times, model_times, d1_mass); d1_mask = np.ma.masked_values(d1_interp, 0)
d2_interp = np.interp(interp_times, model_times, d2_mass); d2_mask = np.ma.masked_values(d2_interp, 0)
d4_interp = np.interp(interp_times, model_times, d4_mass); d4_mask = np.ma.masked_values(d4_interp, 0)

idx1 = np.where(interp_times == 60*60)[0][0]
dm_d1 = d1_mask[idx1] - d1_mask[0] # kg
dm_d2 = d2_mask[idx1] - d2_mask[0]
dm_d4 = d4_mask[idx1] - d4_mask[0]

dias = [dm_d1, dm_d2, dm_d4]

dm_n100 = n100_mass[2] - n100_mass[0] / 1000 # kg/s
dm_n500 = n500_mass[2] - n500_mass[0] / 1000
dm_n1000 = n1000_mass[2] - n1000_mass[0] / 1000

nums = [dm_n100, dm_n500, dm_n1000]
eq_times = np.empty((3,3))

for idia, inum in zip(dias, nums):
    print(idia)
    print(inum)
    eq_times[dias.index(idia), nums.index(inum)] = idia / inum

# %% Add first value to arrays
model_times = np.insert(output_n100['calc_time'], 0, 0)
model_times = model_times/60

temp_dict = output_n100
(n_scwater * moisture_calculations.dropVolume(d_scwater * 1e-3) * constants.rhoI * 1e3) / (spacing ** 3)

n100_mass = output_n100['drop_mass_kg']; n100_mass = np.insert(n100_mass, 0, starting_drop_mass); n100_mass = np.ma.masked_values(n100_mass,0) # n = 100/cm3
n500_mass = output_n500['drop_mass_kg']; n500_mass = np.insert(n500_mass, 0, starting_drop_mass); n500_mass = np.ma.masked_values(n500_mass,0) # n = 500/cm3
n1000_mass = output_n1000['drop_mass_kg']; n1000_mass = np.insert(n1000_mass, 0, starting_drop_mass); n1000_mass = np.ma.masked_values(n1000_mass,0) # n = 500/cm3

r102_mass = output_r102['drop_mass_kg']; r102_mass = np.insert(r102_mass, 0, starting_drop_mass); r102_mass = np.ma.masked_values(r102_mass,0) # n = 100/cm3
r105_mass = output_r105['drop_mass_kg']; r105_mass = np.insert(r105_mass, 0, starting_drop_mass); r105_mass = np.ma.masked_values(r105_mass,0) # n = 100/cm3
r110_mass = output_r110['drop_mass_kg']; r110_mass = np.insert(r110_mass, 0, starting_drop_mass); r110_mass = np.ma.masked_values(r110_mass,0) # n = 100/cm3
r120_mass = output_r120['drop_mass_kg']; r120_mass = np.insert(r120_mass, 0, starting_drop_mass); r120_mass = np.ma.masked_values(r120_mass,0) # n = 100/cm3
r130_mass = output_r130['drop_mass_kg']; r130_mass = np.insert(r130_mass, 0, starting_drop_mass); r130_mass = np.ma.masked_values(r130_mass,0) # n = 100/cm3

# %%
plot_type = 'radius' #'iwc', 'radius'

if plot_type == 'mass':
    ylab = 'Mass [kg]'
elif plot_type == 'iwc':
    ylab = 'IWC [g m$^{-3}$]'
elif plot_type == 'radius':
    ylab = 'Radius [mm]'

xunits = 'mins'
#model_times = model_times / 60

if xunits == 'mins':
    xlab = 'Time [minutes]'
elif xunits == 'seconds':
    xlab = 'Time [seconds]'

# %% Plot model drop mass (deposition)
fig = plt.figure(figsize=(6,6))
ax1 = plt.gca()
if plot_type == 'mass':
    ax1.plot(model_times, r102_mass, label='RH$_{ice}$ = 102.5%')
    ax1.plot(model_times, r105_mass, label='RH$_{ice}$ = 105%')
    ax1.plot(model_times, r110_mass, label='RH$_{ice}$ = 110%')
    ax1.plot(model_times, r120_mass, label='RH$_{ice}$ = 120%')
    ax1.plot(model_times, r130_mass, label='RH$_{ice}$ = 130%')
elif plot_type == 'iwc':
    ax1.plot(model_times, r102_mass*1000/(1000e3), label='RH$_{ice}$ = 102.5%')
    ax1.plot(model_times, r105_mass*1000/(1000e3), label='RH$_{ice}$ = 105%')
    ax1.plot(model_times, r110_mass*1000/(1000e3), label='RH$_{ice}$ = 110%')
    ax1.plot(model_times, r120_mass*1000/(1000e3), label='RH$_{ice}$ = 120%')
    ax1.plot(model_times, r130_mass*1000/(1000e3), label='RH$_{ice}$ = 130%')
elif plot_type == 'radius':
    ax1.plot(model_times, (moisture_calculations.dropDiameter(r102_mass / constants.rhoI) / 1e-3)/2, label='RH$_{ice}$ = 102.5%', color='#01d298', linewidth=3)
    ax1.plot(model_times, (moisture_calculations.dropDiameter(r105_mass / constants.rhoI) / 1e-3)/2, label='RH$_{ice}$ = 105%', color='#04a779', linewidth=3)
    ax1.plot(model_times, (moisture_calculations.dropDiameter(r110_mass / constants.rhoI) / 1e-3)/2, label='RH$_{ice}$ = 110%', color='#087e5c', linewidth=3)
    ax1.plot(model_times, (moisture_calculations.dropDiameter(r120_mass / constants.rhoI) / 1e-3)/2, label='RH$_{ice}$ = 120%', color='#095740', linewidth=3)
    ax1.plot(model_times, (moisture_calculations.dropDiameter(r130_mass / constants.rhoI) / 1e-3)/2, label='RH$_{ice}$ = 130%', color='#073325', linewidth=3)

ax1.set_ylabel(ylab)
ax1.set_xlabel(xlab)
ax1.grid()
ax1.set_title('Vapor Deposition sensitivity', loc='right')
ax1.legend()
plt.show()
#plt.savefig("Q:\\My Drive\\phd\\winter_storms\\seeder_feeder\\symposium\\deposition_sensitivity_{0}.png".format(plot_type), dpi=300, bbox_inches='tight')

# %% Plot model drop mass (riming)
fig = plt.figure(figsize=(6,6))
ax1 = plt.gca()
if plot_type == 'mass':
    ax1.plot(model_times, n100_mass, label='n = 100 cm$^{-3}$', color='#bae4bc')
    ax1.plot(model_times, n500_mass, label='n = 500 cm$^{-3}$', color='#7bccc4')
    ax1.plot(model_times, n1000_mass, label='n = 1000 cm$^{-3}$', color='#2b8cbe')
elif plot_type == 'iwc':
    ax1.plot(model_times, n100_mass*1000/(1000e3), label='n = 100 cm$^{-3}$')
    ax1.plot(model_times, n500_mass*1000/(1000e3), label='n = 500 cm$^{-3}$')
    ax1.plot(model_times, n1000_mass*1000/(1000e3), label='n = 1000 cm$^{-3}$')
elif plot_type == 'radius':
    ax1.plot(model_times, (moisture_calculations.dropDiameter(n100_mass / constants.rhoI) / 1e-3)/2, label='n = 100 cm$^{-3}$', color='#00c6e9', linewidth=3)
    ax1.plot(model_times, (moisture_calculations.dropDiameter(n500_mass / constants.rhoI) / 1e-3)/2, label='n = 500 cm$^{-3}$', color='#008ec3', linewidth=3)
    ax1.plot(model_times, (moisture_calculations.dropDiameter(n1000_mass / constants.rhoI) / 1e-3)/2, label='n = 1000 cm$^{-3}$', color='#075891', linewidth=3)
ax1.set_ylabel(ylab)
ax1.set_xlabel(xlab)
ax1.grid()
ax1.set_title('Riming sensitivity', loc='right')
ax1.legend()
#plt.show()
plt.savefig("Q:\\My Drive\\phd\\winter_storms\\seeder_feeder\\symposium\\riming_sensitivity_{0}.png".format(plot_type), dpi=300, bbox_inches='tight')

#%% both
fig = plt.figure(figsize=(6,6))
ax1 = plt.gca()
ax1.plot(model_times, (moisture_calculations.dropDiameter(n100_mass / constants.rhoI) / 1e-3)/2, label='n = 100 cm$^{-3}$', color='#00c6e9', linewidth=3)
ax1.plot(model_times, (moisture_calculations.dropDiameter(n500_mass / constants.rhoI) / 1e-3)/2, label='n = 500 cm$^{-3}$', color='#008ec3', linewidth=3)
ax1.plot(model_times, (moisture_calculations.dropDiameter(n1000_mass / constants.rhoI) / 1e-3)/2, label='n = 1000 cm$^{-3}$', color='#075891', linewidth=3)
ax1.plot(model_times, (moisture_calculations.dropDiameter(r102_mass / constants.rhoI) / 1e-3) / 2,
         label='RH$_{ice}$ = 102.5%', color='#01d298', linewidth=2)
ax1.plot(model_times, (moisture_calculations.dropDiameter(r105_mass / constants.rhoI) / 1e-3) / 2,
         label='RH$_{ice}$ = 105%', color='#04a779', linewidth=2)
ax1.plot(model_times, (moisture_calculations.dropDiameter(r110_mass / constants.rhoI) / 1e-3) / 2,
         label='RH$_{ice}$ = 110%', color='#087e5c', linewidth=2)
ax1.plot(model_times, (moisture_calculations.dropDiameter(r120_mass / constants.rhoI) / 1e-3) / 2,
         label='RH$_{ice}$ = 120%', color='#095740', linewidth=2)
ax1.plot(model_times, (moisture_calculations.dropDiameter(r130_mass / constants.rhoI) / 1e-3) / 2,
         label='RH$_{ice}$ = 130%', color='#073325', linewidth=2)

ax1.set_ylabel(ylab)
ax1.set_xlabel(xlab)
#ax1.set_ylim([0,60])
ax1.set_yscale('log')
ax1.set_xscale('log')
ax1.grid()
#ax1.set_title('Riming sensitivity', loc='right')
ax1.legend()
plt.show()
#plt.savefig("Q:\\My Drive\\phd\\winter_storms\\seeder_feeder\\symposium\\both_sensitivity_{0}.png".format(plot_type), dpi=300, bbox_inches='tight')

# %% Summary plot
fig = plt.figure(figsize=(8,8))
gs = fig.add_gridspec(ncols=4,nrows=3, wspace=0.1, hspace=0.4)

ax1 = fig.add_subplot(gs[:2, 0]) # RH
ax1.plot(RHi, grid_mids, label='RHice')
ax1.plot(RHw, grid_mids, label='RHwater')
ax1.set_ylabel('Height [m]')
ax1.set_xlabel('Relative Humidity [%]')
ax1.set_ylim([bottom, top])
ax1.set_xlim([60,110])
ax1.grid()
ax1.legend()

ax2 = fig.add_subplot(gs[:2, 1]) # Temperature
ax2.plot(temp, grid_mids, color='red', label='Temperature')
ax2.plot(dewp, grid_mids, color='green', label='Dewpoint')
#ax2.set_ylabel('Height [m]')
ax2.set_xlabel('Temperature [degC]')
ax2.set_ylim([bottom, top])
ax2.set_xlim([-40,10])
ax2.axes.yaxis.set_ticklabels([])
ax2.grid()
ax2.legend()

ax3 = fig.add_subplot(gs[:2, 2]) # Supercooled LWC
ax3.plot((n_scwater * moisture_calculations.dropVolume(d_scwater * 1e-3) * (constants.rhoI/1000))/(1000**3), \
         grid_mids, color='dodgerblue', label='Supercooled \nLWC') # g/m3
#ax3.set_ylabel('Height [m]')
ax3.set_xlabel('LWC [g m$^{-3}$]')
ax3.set_ylim([bottom, top])
ax3.axes.yaxis.set_ticklabels([])
ax3.grid()
ax3.legend()

ax4 = fig.add_subplot(gs[2,:4]) # dm/dt
ax4.plot(model_times, col_mass_flux_masked)#np.ediff1d(col_drop_mass, to_end=0)/time_step)
ax4.set_ylabel('Mass Flux [kg]')
ax4.set_xlabel('Time [s]')
ax4.set_xlim([model_times[0], model_times[-1]])
ax4.grid()

ax5 = fig.add_subplot(gs[:2,3]) # dm/dz
ax5.plot(col_mass_flux_masked, grid_mids[::-1])
#ax5.set_ylabel('Height [m]')
ax5.set_xlabel('Mass Flux [kg]')
ax5.set_ylim([bottom, top])
ax5.axes.yaxis.set_ticklabels([])
ax5.grid()

plt.show()

# %% Plot model drop diameter
col_drop_dia = np.ma.masked_values(col_drop_dia, 0)

fig = plt.figure(figsize=(8, 8))
ax1 = plt.gca()
ax1.plot(col_drop_dia, drop_height)
ax1.set_ylabel('Height [m]')
ax1.set_ylim([0,10000])
ax1.set_xlabel('Diameter [mm]')
ax1.yaxis.set_minor_locator(FixedLocator(grid_tops))
ax1.set_title('RHice = {0}%'.format(str(RHi[0])), loc='right')
ax1.grid()
plt.show()

# %% Plot model drop mass
col_drop_mass = np.ma.masked_values(col_drop_mass, 0)

fig = plt.figure(figsize=(8, 8))
ax1 = plt.gca()
ax1.plot(col_drop_mass, drop_height)
ax1.set_ylabel('Height [m]')
ax1.set_ylim([0,10000])
ax1.set_xlabel('Mass [kg]')
ax1.yaxis.set_minor_locator(FixedLocator(grid_tops))
ax1.grid()
ax1.set_title('RHice = {0}%'.format(str(RHi[0])), loc='right')
plt.show()

# %% Plot model mass_flux
fig = plt.figure(figsize=(8, 8))
ax1 = plt.gca()
ax1.plot(col_mass_flux, grid_mids[::-1])
ax1.set_ylabel('Height [m]')
#ax1.set_xlabel('Diameter [micrometers]')
ax2 = ax1.twinx();
ax2.invert_yaxis()
ax2.plot(col_mass_flux, model_times)
ax2.set_ylabel('Time [s]')
ax1.grid()
plt.show()

# %% Plotting temperature and RHice profiles
fig = plt.figure(figsize=(4, 6))
ax1 = plt.gca()
ax1.plot(temp, grid_mids)
ax1.set_ylabel('Height [m]')
ax1.set_xlabel('Temperature [degC]')
ax2 = ax1.twiny()
ax2.plot(RHi, grid_mids)
ax2.set_xlabel('RHice')
plt.show()

# %% Temp and Dewp plots
fig = plt.figure(figsize=(4, 6))
ax = plt.gca()
ax.plot(temp, grid_mids, color='red', label='Temperature')
ax.plot(dewp, grid_mids, color='green', label='Dewpoint')
ax.set_ylabel('Height [m]')
ax.set_xlabel('Temperature [degC]')
ax.grid()
plt.legend()
plt.show()

# %% RH and RHice plots
fig = plt.figure(figsize=(4, 6))
ax = plt.gca()
ax.plot(RHi, grid_mids, label='RHice')
ax.plot(RHw, grid_mids, label='RHwater')
ax.set_ylabel('Height [m]')
ax.set_xlabel('Relative Humidity [%]')
ax.grid()
plt.legend()
plt.show()

# %% Fall speed plots
fig = plt.figure(figsize=(4, 6))
ax = plt.gca()
ax.plot(moisture_calculations.fallspeedCorrection(1, rho_mids), grid_mids)
ax.set_ylabel('Height [m]')
ax.set_xlabel('Fall speed [m/s]')
ax.set_ylim([0,10000])
ax.grid()
plt.show()

# %%
temps = np.array([-40, -30, -20, -10, 0, 10, 20, 30])
Ka = np.array([2.07e-2, 2.16e-2, 2.24e-2, 2.32e-2, 2.4e-2, 2.48e-2, 2.55e-2, 2.63e-2])
mk, bk = np.polyfit(temps, Ka, 1)
Dv = np.array([1.62e-5, 1.76e-5, 1.91e-5, 2.06e-5, 2.21e-5, 2.36e-5, 2.52e-5, 2.69e-5])
md, bd = np.polyfit(temps, Dv, 1)
plt.figure()
plt.scatter(temps, Ka)
plt.plot(temps, (mk * temps + bk))
plt.show()

plt.figure()
plt.scatter(temps, Dv)
plt.plot(temps, (md * temps + bd))
plt.show()
#%% # set starting values
iheight = grid['top']
idiameter = starting_dsd['diameter']
idropmass = starting_dsd['mass']

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
        break

    height_idx = moisture_calculations.nearest_ind(grid['grid_mids'],
                                                   iheight)  # np.where(grid_tops == iheight)[0][0]

    # Find supersaturation, Fk, and Fd at height
    supersat = (environment['e'][height_idx] / environment['esi'][height_idx]) - 1
    Fk = moisture_calculations.FkCalc(environment['temp_C'][height_idx], iceFlag=True)  # [m s kg-1]
    Fd = moisture_calculations.FdCalc(environment['temp_C'][height_idx], environment['esi'][height_idx])  # [m s kg-1]

    # Calculate fall speed
    fall_speed = model_calculations.fall_speed_from_d(idiameter)
    fall_speed_corrected = model_calculations.fallspeedCorrection(fall_speed,
                                                                  environment['air_density'][height_idx])

    # Calculate sublimation mass flux
    sublimation_flux = (4 * np.pi * (idiameter / 2) * 1e-3 * supersat) / (Fk + Fd)  # dm/dt [kg s-1]

    # Calculate riming mass flux
    A = np.pi * (((idiameter / 2) + (environment['scwater']['d'][height_idx] / 2)) * 1e-3) ** 2
    riming_flux = A * abs(fall_speed_corrected - 0) * environment['scwater']['coll_eff'] * \
                  (environment['scwater']['n'][height_idx] * moisture_calculations.dropVolume(
                      environment['scwater']['d'][height_idx] * 1e-3) * constants.rhoI) / (
                          grid['spacing'] * 1000 ** 2)

    #print(str(riming_flux))

    # update particle
    total_flux = (sublimation_flux + riming_flux) * time_step
    new_dropmass = idropmass + total_flux  # [kg]
    if new_dropmass <= 0:
        print("Drop mass <= 0; exiting loop")
        break
    new_diameter = moisture_calculations.dropDiameter(idropmass / constants.rhoI) / 1e-3  # [mm]

    # store variables
    col_mass_flux[time_idx] = total_flux
    col_drop_mass[time_idx] = new_dropmass
    col_drop_dia[time_idx] = new_diameter
    col_latent_e[time_idx] = sublimation_flux * time_step * constants.Ls
    col_iwc[time_idx] = new_dropmass / 1000e3
    drop_height[time_idx] = iheight

    # update for next iteration
    idropmass = new_dropmass
    idiameter = new_diameter
    iheight = iheight - time_step * fall_speed_corrected

    if iheight <= 0:
        print("Height <= 0; exiting loop")
        break

output_r102 = {
    'calc_time': model_times,
    'mass_flux': col_mass_flux,
    'drop_mass_kg': col_drop_mass,
    'drop_dia': col_drop_dia,
    'latent_e': col_latent_e,
    'iwc': col_iwc}

