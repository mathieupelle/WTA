## Unsteady BEMT code (Dynamic inflow)
## Assignment A - Wind Turbine Aeroelasticity
## TU Delft 2021


#%% Input necessary packages 
import numpy as np
import matplotlib.pyplot as plt
from utilities_uBEMT import Rotor,BEMT

#%% Case A: Dynamic inflow due to change in rotor configuration  
# A.1 Step change in thrust coefficient

Geometry = Rotor(N_radial_sections = 30) #Define the rotor geometry

#Define which steps do we want to analyse. This will be the summary variable for all the case
CT_step = {'cases':[[0.5,0.9],
           [0.9,0.5],
           [0.2,1.1],
           [1.1,0.4]],
           'pitch':[],
           'results':[]}

Calc = BEMT(Geometry) #Firstly initialize the BEMT class to compute the calculations
Calc.CpLambda(TSR_list = [10], theta_list = list(np.linspace(-7,5))) #Calculate the Cp/Ct-theta-tsr contours

#Loop through each of these cases
for i,val in enumerate(CT_step['cases']):
    #Get the pitch angle for each value of the step
    pitch_angle = [Calc.getPitchAngle_fromCT(CT = val[0], TSR = 10), Calc.getPitchAngle_fromCT(CT = val[1], TSR = 10)]
    
    #Store the pitch values in the dictionary of the CT_step summary
    CT_step['pitch'].append(pitch_angle)
    
    #Build conditions dictionary necessary for the calculation of the unsteady BEMT
    time_arr = np.linspace(0,10)
    cond = {'wind_speed': 10*np.ones(len(time_arr)),
        'pitch_angle': np.concatenate((np.array(pitch_angle),pitch_angle[1]*np.ones(len(time_arr)-2)),axis=None),
        'yaw_angle': np.zeros(len(time_arr))}
    
    #Run BEMT
    Calc.Solver(time = time_arr, conditions = cond, DI_Model = "PP")
    
    #Store the results in the summary dictionary
    CT_step['results'].append(Calc.Results)
    
    
#%% A.2 - Sinusoidal change in quasi-steady thrust coefficient

#Define time array and reduced frequency
time_arr = np.arange(0,10,0.01)
omega = np.arange(0.05,0.301,0.05)






