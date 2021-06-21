"""
Created on Mon Jun 21 21:27:41 2021

@author: nilsg

For testing Øye dynamic inflow method BEMT.Oye() from utilities_DI.py
"""

import numpy as np
from utilities_DI import Rotor,BEMT

Test = Rotor(N_radial_sections = 40)
Test.SetOperationalData(wind_speed = 10, TSR = 10, yaw = 0)
Sol = BEMT(Test)

# we will now calculate and plot the solution of a step change in thrust coefficient using the Oye model

# we define the value of U_infinity and the radius of the actuator
Uinf=1 # U_infinity
R=1 # radius of the actuator

# define time array
dt=0.005 # we define the time step
time=np.arange(0, 20, dt) # we create the array "time"


# define Ct and induction at  t0
Ct0=np.array([.5])  # this it the value of thrust coefficient at time t0, the start of the time array
vind0=(Ct0)*Uinf # this is the quasi-steady value of induction at time t0, calculated from Ct0

# define quasi-steady solution of Ct and induction at  t>=t1
Ct1=np.array([-0.85]) # this it the value of thrust coefficient at time t1
vind1=-ainduction(-Ct1)*Uinf # this is the quasi-steady value of induction at time t1, calculated from Ct1

# define Ct as a function of time
Ct = np.zeros(np.shape(time))+Ct0 # we initialize the array of thrust coefficient, setting all initial values at Ct0

# change Ct for time above t1
t1=5 # we define t1, when Ct experiences a step change
Ct[time>=t1] = Ct1 # we correct the values of Ct for the instants that time is after t1, to a value of Ct1. We therefore 
                    # define the step change from Ct0 to Ct1
    
#set arrays for induction
vind = np.zeros(np.shape(time)) # we create the array of induction velocity
vind[0]=vind0 # we initialize the first value to the induction velocity at t0, calculated from Ct0

# the Oye model requires the use of a temporary value of induction vint. Here, we initialize the value of vint
vint=vind0


# solve the equation in time of the value of induction by using the Oye model 
for i,val in enumerate(time[:-1]):
    vind[i+1],vint=oye_dynamic_inflow(vind[i], Ct[i], Ct[i+1], vint, Uinf, R, 0.95*R,dt)  