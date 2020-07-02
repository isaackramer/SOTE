"""SOTE script"""
import numpy as np
import SOTE_class
import math
import pandas as pd

np.random.seed(10) # set seed

years = 2 # simulation length in years
dt = 0.05 # time step in days
time = np.arange(0, years*365+dt, dt)

# instantiate class and function (soil_type must be one of: 'class_1', 'class_2', 'class_3')
eq = SOTE_class.SOTE(Cirr=10, Eirr=0.55, C_init=20.0, E_init=0.02, ET_w = .1,
                     s_init = 0.3, depth = 350, soil_type = 'class_1', Crain = 0.0,
                     Erain = 0.0, rain_prob = 0.3, mean_height = 10, days = 130,
                     dt=dt, ET_ratio = 1.1, runs = 100)

# bounds of rainy season
rain_start = 365 - eq.par['days']
rain_end = 365
rain_mid = 365 - eq.par['days']/2
sin_shift = rain_mid - 273.5

# inital values based on class initalization
s0 = eq.par['s_init']          # soil moisture [nondimensional]
w0 = s0 * eq.soil_parms['nZr'] # water content [mm]
C0 = eq.par['C_init']          # salt concentration [mmol_c/L]
q0 = eq.par['C_init']*w0       # salt mass content [mmol_c]
E0 = eq.par['E_init']          # exchangeable sodium fraction [nondimensional]
Krel0 = 1                      # relative Ks [nondimensional]

# we will run an ensemble of N=runs simulations at once
# vectorize inital values
s = np.full((eq.par['runs'],), s0)
w = np.full((eq.par['runs'],), w0)
q = np.full((eq.par['runs'],), q0)
C = np.full((eq.par['runs'],), C0)
E = np.full((eq.par['runs'],), E0)
Krel = np.full((eq.par['runs'],), Krel0)

# array to record average results (t, s, C, E, Krel)
results_avg = np.array([[time[0], s0, C0, E0, Krel0]])

# main loop
for tim in time[1:]:
    # ET_max dependent on time of year
    eq.ETmax = 2 * np.sin((tim - sin_shift) * ((2 * np.pi) / 365)) + 5

    # stochastic rainfall contribution
    DOY = (math.floor(tim) % 365) # day of year since 1 January
    if (DOY >= rain_start and DOY <= rain_end):
        # rainy season
        eq.rain_height()
        eq.Irr = 0.0
    else:
        # dry season
        eq.event_height = np.zeros(eq.par['runs'],)
        eq.rain_rate = np.zeros(eq.par['runs'],)
        eq.Irr = eq.par['ET_ratio'] * eq.ETmax

    # calculates net change in water content + salinity and sodicity of input and output water
    eq.water_net(s, q/w, E)

    # integration
    eq.derivs()
    new_values = eq.rk4_step([q, E, s])
    q, E, s = new_values[0,:]
    # water content must be less than 1 (saturation excess)
    s = np.where(s > 1.0, 1.0, s)
    w = s*eq.soil_parms['nZr']

    # update array to track results (t, s, C, E, Krel)
    step = np.array([[tim,
                      np.mean(s),
                      np.mean(q)/np.mean(w),
                      np.mean(E),
                      np.mean(eq.Ksat/eq.soil_parms['Ks'])]])
    results_avg = np.concatenate((results_avg, step))

# save output
results_avg = pd.DataFrame(results_avg,
                   columns = ['time (days)',
                              'relative soil water content (nondim)',
                              'electrolyte concentration (mmol_c/L)',
                              'exchangeable sodium fraction (nondim)',
                              'relative Ks (nondim)'])
results_avg.to_csv('results_average.csv', index=False)

ends = np.array([s, q/w, E, eq.Ksat/eq.soil_parms['Ks']])
np.savetxt('final_values.csv', ends, delimiter=",")
