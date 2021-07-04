## Unsteady BEMT code (Dynamic inflow)
## Assignment A - Wind Turbine Aeroelasticity
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

 
#Define the models that we want to analyse
DI_models = ['Steady','PP','LM','O']

Calc = BEMT(Geometry) #Initialize the BEMT class to compute the calculations
Calc.CpLambda(TSR_list = list(np.linspace(5,15,50)), theta_list = list(np.linspace(-15,10,100))) #Calculate the Cp/Ct-theta-tsr contours

#Plot contours
PlotContours(Calc.Cp_lambda)


#%% Case A: Dynamic inflow due to change in rotor configuration  
# A.1 Step change in thrust coefficient

#Define which steps do we want to analyse. This will be the summary variable for all the case
CT_step = {'cases':[[0.5,0.9],
           [0.9,0.5],
           [0.2,1.1],
           [1.1,0.4]],
           'pitch':[]}

#Initialise the results variable
CT_step['results'] = [[],[],[],[]]

# Loop through each of these cases
for i,val in enumerate(CT_step['cases']):
    #Get the pitch angle for each value of the step
    pitch_angle = [Calc.getPitchAngle_fromCT(CT = val[0], TSR = 10), Calc.getPitchAngle_fromCT(CT = val[1], TSR = 10)]
    
    #Store the pitch values in the dictionary of the CT_step summary
    CT_step['pitch'].append(pitch_angle)
    
    #Build conditions dictionary necessary for the calculation of the unsteady BEMT

    time_arr = np.arange(0,30.05,0.05)
    CT_step['time'] = time_arr
    cond = {'wind_speed': 10*np.ones(len(time_arr)),
        'pitch_angle': np.concatenate((np.array(pitch_angle),pitch_angle[1]*np.ones(len(time_arr)-2)),axis=None),
        'yaw_angle': np.zeros(len(time_arr)),
        'TSR': 10*np.ones(len(time_arr))}
    
    #cond['pitch_angle'] = 0*np.ones(len(time_arr))
    
    #Run BEMT for each model
    for j,model in enumerate(DI_models):
        print('Running case',i*len(DI_models)+j+1,'out of',len(CT_step['cases'])*len(DI_models))

        Calc.Solver(time = time_arr, conditions = cond, DI_Model = model)
    
        #Store the results in the summary dictionary
        CT_step['results'][j].append(Calc.Results)
    
#Save the results
if saving:
    file = open("CT_step.pkl","wb")
    pickle.dump(CT_step,file)
    file.close()    
    
#%% A.2 - Sinusoidal change in quasi-steady thrust coefficient

# Define time array and reduced frequency
time_arr = np.arange(0,130,0.2)
omega = [0.05,0.15,0.3]

# Input the sinusoisdal scenarios that will be explored
CT_sin = {'CT_0': [.5,.9,.2],
          'Delta_CT': [.5,.3,.7],
          'pitch_arr': np.zeros((3,len(time_arr))),
          'CT_arr': np.zeros((3,len(time_arr))),
          'results': [[[None]*3]*len(omega)]*4,
          'omega': omega}

#Create the empty results list
CT_sin['results'] = [[[[],[],[]] for guillem in range(len(omega))] for mathieu in range(len(DI_models))]

#Store time array
CT_sin['time'] = time_arr


# Explore all the frequencies defined in omega
for i,val in enumerate(omega):
    #Within each frequency, we want to evaluate all the cases listed in CT_sin
    for j in range(len(CT_sin['CT_0'])):
        #Start by defining the thrust sinusoidal wave
        CT = CT_sin['CT_0'][j] + CT_sin['Delta_CT'][j]*np.sin(val*time_arr)
        #Convert it to pitch 
        pitch_arr = Calc.getPitchAngle_fromCT(CT = CT,TSR = 10)
        #Store it in the summary dictionary
        CT_sin['pitch_arr'][j,:] = pitch_arr
        CT_sin['CT_arr'][j,:] = CT
        
        #Build the conditions dict necessary for the unsteady BEMT
        cond = {'wind_speed': 10*np.ones(len(time_arr)),
                'pitch_angle': pitch_arr,
                'yaw_angle': np.zeros(len(time_arr)),
                'TSR': 10*np.ones(len(time_arr))}
        
        #Run BEMT for each DI model
        for k,model in enumerate(DI_models):
            print('Running case',i*(len(CT_sin['CT_0'])*len(DI_models))+j*(len(DI_models))+k+1,'out of',len(CT_sin['CT_0'])*len(omega)*len(DI_models))
            Calc.Solver(time = time_arr, conditions = cond, DI_Model = model)
            
            #Store the results
            CT_sin['results'][k][i][j] = Calc.Results
                  

if saving:
    file = open("CT_sin.pkl","wb")
    pickle.dump(CT_sin,file)
    file.close()    


#%% B.1 Dynamic inflow due to change in wind speed - Step change in Uinf (keep the same rotational speed)

## INPUTS
#Baseline velocity and tip speed ratio for the calculations
U_0 = 10 
TSR_0 = 10

#Velocity steps
U_step = {'cases':[[1,1.5],
           [1,0.7],
           [1,1.2],
           [1,0.9]]}

#Calculate baseline rotor speed
Omega_0 = TSR_0*U_0/Geometry.radius

#Calculate the pitch angle corresponding to the optimal CT=8/9
pitch =  Calc.getPitchAngle_fromCT(CT = 8/9,TSR = TSR_0)

#Create time array
time_arr = np.linspace(0,1,20)
U_step['time'] = time_arr

#Create results array
U_step['results'] = [[[],[],[],[]] for nils in range(len(DI_models))]


## Iterate for each case
for i,step in enumerate(U_step['cases']):
    #Firstly, calculate the new TSR value
    TSR_step = Omega_0*Geometry.radius/(U_0*step[1])
    
    #Build conditions dictionary
    cond = {'wind_speed': np.concatenate((np.array(U_0),U_0*step[1]*np.ones(len(time_arr)-1)),axis=None),
        'pitch_angle': pitch*np.ones(len(time_arr)),
        'yaw_angle': np.zeros(len(time_arr)),
        'TSR': np.concatenate((np.array(TSR_0),TSR_step*np.ones(len(time_arr)-1)),axis=None)}   
    
    #Run BEMT for each model
    for j,model in enumerate(DI_models):
        print('U-step: running case',i*len(DI_models)+j+1,'out of',len(U_step['cases'])*len(DI_models))
        Calc.Solver(time = time_arr, conditions = cond, DI_Model = model) 
        
        #Store the results
        U_step['results'][j][i] = Calc.Results        
 
if saving:
    file = open("U_step.pkl","wb")
    pickle.dump(U_step,file)
    file.close()

#%% B.2 Dynamic inflow due to change in wind speed - Sinusoidal space

## INPUTS
#Baseline velocity and tip speed ratio for the calculations
U_0 = 10 
TSR_0 = 10

#Velocity steps
U_sin = {'U_1':[1,0.7,1.2],
         'Delta_U':[.5,.3,.5],
         'omega':np.arange(0.05,0.301,0.05)}

#Create time array
time_arr = np.arange(0,10,0.1)
U_sin['time'] = time_arr

#Create empty arrays for U and TSR
U_sin['U_arr'] = np.zeros((len(U_sin['U_1']),len(time_arr)))
U_sin['TSR_arr'] = np.zeros((len(U_sin['U_1']),len(time_arr)))

#Calculate baseline rotor speed
Omega_0 = TSR_0*U_0/Geometry.radius

#Calculate the pitch angle corresponding to the optimal CT=8/9
pitch =  Calc.getPitchAngle_fromCT(CT = 8/9,TSR = TSR_0)

#Create results array
U_sin['results'] = [[[[],[],[]] for guillem in range(len(omega))] for mathieu in range(len(DI_models))]

## Iterate for each case
for i,freq in enumerate(omega):
    for j in range(len(U_sin['U_1'])):
        #Start by defining the velocity sinusoidal wave
        U_arr = (U_sin['U_1'][j] + U_sin['Delta_U'][j]*np.sin(freq*time_arr))*U_0
        #Convert it to tip speed ratio 
        TSR_arr = Omega_0*Geometry.radius/U_arr
        #Store it in the summary dictionary
        U_sin['U_arr'][j,:] = U_arr
        U_sin['TSR_arr'][j,:] = TSR_arr
        
        #Build the conditions dict necessary for the unsteady BEMT
        cond = {'wind_speed': U_arr,
                'pitch_angle': pitch*np.ones(len(time_arr)),
                'yaw_angle': np.zeros(len(time_arr)),
                'TSR': TSR_arr}
        
        #Run BEMT for each DI model
        for k,model in enumerate(DI_models):
            print('Running case',i*(len(U_sin['U_1'])*len(DI_models))+j*(len(DI_models))+k+1,'out of',len(U_sin['U_1'])*len(omega)*len(DI_models))
            Calc.Solver(time = time_arr, conditions = cond, DI_Model = model)
            
            #Store the results
            U_sin['results'][k][i][j] = Calc.Results
            
if saving:
    file = open("U_sin.pkl","wb")
    pickle.dump(U_sin,file)
    file.close()     
        
#%% Studying U step with the same cases as the CT step

#Build the cases
U_step_comp = {'cases':[.5,.9],
               'pitch':[],
               'TSR':[8,np.nan],
               'results':[[],[],[],[]]}

#Set the baseline wind speed and TSR
U_0 = 10 
TSR_0 = U_step_comp['TSR'][0]

#Calculate baseline rotor speed
Omega_0 = TSR_0*U_0/Geometry.radius

#Create time array
time_arr = np.arange(0,30.1,0.1)
U_step_comp['time'] = time_arr

#Calculate the pitch angle corresponding to the CT of case A.1
U_step_comp['pitch'] = Calc.getPitchAngle_fromCT(CT = U_step_comp['cases'][0], TSR = TSR_0)

#Calculate the necessary TSR to go to the next step
U_step_comp['TSR'][1] = Calc.getTSR_fromCT(CT = U_step_comp['cases'][1],theta = U_step_comp['pitch'])

#Calculate U that gives the necessary TSR
U_1 = Omega_0*Geometry.radius/U_step_comp['TSR'][1]

U_step_comp['U'] = [U_0,U_1]

#Build the conditions dictionary
cond = {'wind_speed': np.concatenate((np.array(U_0),U_1*np.ones(len(time_arr)-1)),axis=None),
        'pitch_angle': U_step_comp['pitch']*np.ones(len(time_arr)),
        'yaw_angle': np.zeros(len(time_arr)),
        'TSR': np.concatenate((np.array(TSR_0),U_step_comp['TSR'][1]*np.ones(len(time_arr)-1)),axis=None)}

#Run uBEMT for each 
for i,model in enumerate(DI_models):
        print('Running case',i+1,'out of',len(DI_models))
        Calc.Solver(time = time_arr, conditions = cond, DI_Model = model)
    
        #Store the results in the summary dictionary
        U_step_comp['results'][i] = Calc.Results 

if saving:
    file = open("U_step_comp.pkl","wb")
    pickle.dump(U_step_comp,file)
    file.close()   
