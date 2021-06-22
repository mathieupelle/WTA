"""
Created on Mon Jun 21 21:27:41 2021

@author: nilsg

For testing Øye dynamic inflow method BEMT.Oye() from utilities_DI.py
"""

import numpy as np
from utilities_uBEMT import Rotor,BEMT
import matplotlib.pyplot as plt

# create rotor geometry
Test = Rotor(N_radial_sections = 40)
Sol = BEMT(Test)

#%% we will now calculate and plot the solution of a step change in thrust coefficient using the Oye model

# we define the value of U_infinity and the radius of the actuator
Test.SetOperationalData(wind_speed = 1, TSR = 1, yaw = 0) # wind_speed
Test.radius = 1 # radius of the actuator

# define time array
dt = 0.005 # we define the time step
time = np.arange(0, 20, dt) # we create the array "time"

# define Ct and induction at  t0
Ct0 = np.array([.5])  # this it the value of thrust coefficient at time t0, the start of the time array
a0 = Sol.NewInductionFactor(CT=Ct0) # this is the quasi-steady value of induction factor at time t0, calculated from Ct0
vind0 = Sol.NewInductionFactor(CT=Ct0) * Test.wind_speed # this is the quasi-steady value of induction factor at time t0, calculated from Ct0

# define quasi-steady solution of Ct and induction at  t>=t1
Ct1 = np.array([0.85]) # this it the value of thrust coefficient at time t1
a1 = Sol.NewInductionFactor(CT=Ct1) # this is the quasi-steady value of induction factor at time t1, calculated from Ct1
vind1 = Sol.NewInductionFactor(CT=Ct1) * Test.wind_speed # this is the quasi-steady value of induction velocity at time t0, calculated from Ct0

# define Ct as a function of time
Ct = np.zeros(np.shape(time))+Ct0 # we initialize the array of thrust coefficient, setting all initial values at Ct0

# change Ct for time above t1
t1 = 5 # we define t1, when Ct experiences a step change
Ct[time >= t1] = Ct1 # we correct the values of Ct for the instants that time is after t1, to a value of Ct1. We therefore define the step change from Ct0 to Ct1
    
# set arrays for induction
a = np.zeros(np.shape(time)) # we create the array of induction velocity
a[0] = a0 # we initialize the first value to the induction factor at t0, calculated from Ct0
# vind = np.zeros(np.shape(time)) # we create the array of induction velocity
# vind[0] = vind0 # we initialize the first value to the induction velocity at t0, calculated from Ct0

# the Oye model requires the use of a temporary value of induction vint. Here, we initialize the value of vint
vint = a[0] * Test.wind_speed
# vint = vind[0]

#%% # solve the equation in time of the value of induction by using the Oye model
for i,val in enumerate(time[1:]):
    # Carlos' version
    # vind[i+1],vint=oye_dynamic_inflow(vind[i], Ct[i], Ct[i+1], vint, Uinf, R, 0.95*R,dt)
    # def oye_dynamic_inflow(vz, Ct1, Ct2, vint, Uinf, R, r,dt,glauert=False)
    
    # Our version
    # a[i] = vind[i] / Test.wind_speed
    # vind[i+1], vint = Sol.Oye(a[i], Ct[i], Ct[i+1], vint, Test.radius, dt, Test.wind_speed)
    i += 1
    a[i], vint = Sol.Oye(a[i-1], Ct[i-1], Ct[i], vint, Test.radius, dt, Test.wind_speed)
    # def Oye(a1, CT1, CT2, vint1, r, dt, wind_speed)
    
#%% Plot
fig,ax1 = plt.subplots()
ax2 = ax1.twinx()

lns1 = ax1.plot((time-t1)*Test.radius/Test.wind_speed, (Ct-Ct0)/(Ct1-Ct0), color='green',linestyle='dashdot', label=r'$C_t$')
lns2 = ax2.plot((time-t1)*Test.radius/Test.wind_speed, (a-a0)/(a1-a0), color='blue',linestyle='-', label= r'$a$')
# lns2 = ax2.plot((time-t1)*Test.radius/Test.wind_speed, (vind-vind0)/(vind1-vind0), color='blue',linestyle='-', label= r'$a$')

# define properties of the primary y-axis
# ax1.set_aspect(aspect=20.0) # set aspect ratio of the figure
ax1.set_xlabel(r'$t \frac{R}{U_\infty}$') # label of x-axis
ax1.set_ylabel(r'$\frac{C_T-C_{T_0}}{C_{T_1}-C_{T_0}}$', 
               color='green', fontsize=20) # label of y-axis
ax1.set_xlim(-t1,10) # set limits of x-axis   
ax1.set_ylim(-0.05,1.05) # set limits of x-axis   
ax1.tick_params(axis='y', labelcolor='green') # set the color of the axis

# define properties of secondary axis
ax2.set_ylabel(r'$\frac{a-a_{\left(C_{T_0}\right)}}{a_{\left(C_{T_1}\right)}-a_{\left(C_{T_0}\right)}}$',
               color='blue', fontsize=20) # label of y-axis
ax1.set_ylim(-0.05,1.05) # set limits of x-axis   
ax2.tick_params(axis='y', labelcolor='blue')# set the color of the axis

# here we plot the legend of the figure
lns = lns1+lns2 # add legends
labs = [l.get_label() for l in lns] # extract labels
plt.legend(lns, labs, loc='lower right') # plot legend

ax1.grid(axis='both',which='both') # add a grid

plt.tight_layout() # tight layout to avoid cropping labels
plt.show() # show figure