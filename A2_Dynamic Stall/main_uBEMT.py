## Unsteady BEMT code (Dynamic inflow and dynamic stall)
## Assignment B - Wind Turbine Aeroelasticity
## TU Delft 2021


#%% Input necessary packages
import numpy as np
import matplotlib.pyplot as plt
from utilities_uBEMT import Rotor,BEMT,PlotContours
import pickle

#%% General inputs
#Save results?
saving = True

Geometry = Rotor(N_radial_sections = 15) #Define the rotor geometry

Calc = BEMT(Geometry) #Initialize the BEMT class to compute the calculations
Calc.CpLambda(TSR_list = list(np.linspace(5,15,50)), theta_list = list(np.linspace(-7,5,50))) #Calculate the Cp/Ct-theta-tsr contours

#Plot contours
#PlotContours(Calc.Cp_lambda)

#%% Case Dyn1: Dynamic inflow due to change in inflow velocity azimuthal direction

#Operating conditions
U_0 = 10
TSR_0 = 8
Omega_0 = TSR_0*U_0/Geometry.radius
pitch =  Calc.getPitchAngle_fromCT(CT = 8/9,TSR = TSR_0) #Calculate the pitch angle corresponding to the optimal CT=8/9

#Time and frequency settings
k = 0
omega = k*U_0/Geometry.radius
#T = 2*np.pi/omega
T = 20
time_arr = np.arange(0,5*T,0.1)

#Inflow velocity vector time and and azimuth dependent
N_azimuth = 15
phi = np.linspace(0,2*np.pi,N_azimuth)
U_arr = np.zeros((N_azimuth-1, len(time_arr)))
for i in range(len(phi)-1):
    azimuth = (phi[i]+phi[i+1])/2
    U_arr[i,:] = (1 + 0.5*np.cos(omega*time_arr)*np.cos(azimuth))*U_0

#Models to test
DS_models = ['Steady','BL_noLEsep','BL']

#Initialise the results variable
Dyn1 = {'time':time_arr, 'results': []}

#Build conditions dictionary necessary for the calculation of the unsteady BEMT
cond = {'wind_speed': U_arr,
        'pitch_angle': pitch*np.ones(len(time_arr)),
        'yaw_angle': np.zeros(len(time_arr)),
        'TSR': Omega_0*Geometry.radius/(1 + 0.5*np.cos(omega*time_arr)*U_0),
        'omega': Omega_0*np.ones(len(time_arr))}

#Run BEMT for each model
for i in range(len(DS_models)):
    print('Running case', i+1,'of', len(DS_models))

    #Run simulation
    Calc.Solver(time = time_arr, conditions = cond, DI_Model = 'PP', DS_Model=DS_models[i])

    #Store the results in the summary dictionary
    Dyn1['results'].append(Calc.Results)

#Save the results
if saving:
    file = open("Dyn1.pkl","wb")
    pickle.dump(Dyn1,file)
    file.close()

#%% Case Dyn2: Dynamic inflow due to change in inflow velocity in time and azimuthal direction.

#Operating conditions
U_0 = 10
TSR_0 = 8
Omega_0 = TSR_0*U_0/Geometry.radius
pitch =  Calc.getPitchAngle_fromCT(CT = 8/9,TSR = TSR_0) #Calculate the pitch angle corresponding to the optimal CT=8/9

#Time and frequency settings
k = 0.3
omega = k*U_0/Geometry.radius
T = 2*np.pi/omega
time_arr = np.arange(0,5*T,0.1)

#Inflow velocity vector time and and azimuth dependent
N_azimuth = 15
phi = np.linspace(0,2*np.pi,N_azimuth)
U_arr = np.zeros((N_azimuth-1, len(time_arr)))
for i in range(len(phi)-1):
    azimuth = (phi[i]+phi[i+1])/2
    U_arr[i,:] = (1 + 0.5*np.cos(omega*time_arr)*np.cos(azimuth))*U_0

#Models to test
DS_models = ['BL']
#Initialise the results variable
Dyn2 = {'time':time_arr, 'results': []}

#Build conditions dictionary necessary for the calculation of the unsteady BEMT
cond = {'wind_speed': U_arr,
        'pitch_angle': pitch*np.ones(len(time_arr)),
        'yaw_angle': np.zeros(len(time_arr)),
        'TSR': Omega_0*Geometry.radius/(1 + 0.5*np.cos(omega*time_arr)*U_0),
        'omega': Omega_0*np.ones(len(time_arr))}

#Run BEMT for each model
for i in range(len(DS_models)):
    print('Running case', i+1,'of', len(DS_models))

    #Run simulation
    Calc.Solver(time = time_arr, conditions = cond, DI_Model = 'PP', DS_Model=DS_models[i])

    #Store the results in the summary dictionary
    Dyn2['results'].append(Calc.Results)

#Save the results
if saving:
    file = open("Dyn2.pkl","wb")
    pickle.dump(Dyn2,file)
    file.close()










