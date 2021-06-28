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
CT_step = {'cases':[[0.5,0.5],
           [0.9,0.5],
           [0.2,1.1],
           [1.1,0.4]],
           'pitch':[],
           'results':[]}

CT_step = {'cases':[[0.5,0.9]],
           'pitch':[],
           'results':[]}
Calc = BEMT(Geometry) #Firstly initialize the BEMT class to compute the calculations
Calc.CpLambda(TSR_list = [10], theta_list = list(np.linspace(-7,5))) #Calculate the Cp/Ct-theta-tsr contours

#%% Loop through each of these cases
CT_step['results'] = []
for i,val in enumerate(CT_step['cases']):
    #Get the pitch angle for each value of the step
    pitch_angle = [Calc.getPitchAngle_fromCT(CT = val[0], TSR = 10), Calc.getPitchAngle_fromCT(CT = val[1], TSR = 10)]
    
    #Store the pitch values in the dictionary of the CT_step summary
    CT_step['pitch'].append(pitch_angle)
    
    #Build conditions dictionary necessary for the calculation of the unsteady BEMT
    time_arr = np.linspace(0,1,20)
    cond = {'wind_speed': 10*np.ones(len(time_arr)),
        'pitch_angle': np.concatenate((np.array(pitch_angle),pitch_angle[1]*np.ones(len(time_arr)-2)),axis=None),
        'yaw_angle': np.zeros(len(time_arr))}
    
    #cond['pitch_angle'] = 0*np.ones(len(time_arr))
    
    #Run BEMT
    print('Running case',i+1,'out of',len(CT_step['cases']))
    Calc.Solver(time = time_arr, conditions = cond, DI_Model = "PP")
    
    #Store the results in the summary dictionary
    CT_step['results'].append(Calc.Results)
    
test = CT_step['results'][0]
plt.plot(time_arr,np.mean(test.a,axis=0)[0,:])
plt.plot(time_arr,test.CT)
plt.legend(['$a$','$C_T$'])

plt.plot(test.mu[:,0,0],test.alpha[:,0,:])
    
    
#%% A.2 - Sinusoidal change in quasi-steady thrust coefficient

# Define time array and reduced frequency
time_arr = np.arange(0,0.1,0.01)
omega = np.arange(0.05,0.301,0.05)

# Input the sinusoisdal scenarios that will be explored
CT_sin = {'CT_0': [.5,.9,.2],
          'Delta_CT': [.5,.3,.7],
          'pitch_arr': np.zeros((3,len(time_arr))),
          'CT_arr': np.zeros((3,len(time_arr))),
          'results': [[[None]*3]*len(omega)]*4,
          'omega': omega}

#Select the DI models to test
DI_models = ['Steady','PP','LM','O']

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
                'yaw_angle': np.zeros(len(time_arr))}
        
        #Run BEMT for each DI model
        for k,model in enumerate(DI_models):
            print('Running case',i*(len(CT_sin['CT_0']))+j*(len(DI_models))+k+1,'out of',len(CT_sin['CT_0'])*len(omega)*len(DI_models))
            Calc.Solver(time = time_arr, conditions = cond, DI_Model = model)
            
            #Store the results
            CT_sin['results'][k][i][j] = Calc.Results
            
        







