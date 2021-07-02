# -*- coding: utf-8 -*-
"""
Created on Fri Jul  2 15:17:57 2021

@author: Mathieu Pellé
"""


import numpy as np
import matplotlib.pyplot as plt
import pickle

x = 6
plt.rc('figure', figsize=[46.82 * .5**(.5 * x), 33.11 * .5**(.5 * x)]   )
plt.rc('font', family='serif')

plot_input = {'Cases':['A.1.3', 'A.1.2'], 'Models':['PP'], 'Dimension':['2D'], 'Variable':['CT']}

#Steps in CT
# A.1.1 0.5 --> 0.9
# A.1.2 0.9 --> 0.5
# A.1.3 0.2 --> 1.1
# A.1.4 1.1 --> 0.4

#Oscillations in CT
# A.x.x_omega --> omega = desired frequency
# A.2.1 --> 0.5 to 0.5
# A.2.2 --> 0.9 to 0.3
# A.2.3 --> 0.2 to 0.7


dic = plot_input

cases = dic['Cases']
models =  dic['Models']
dimension = dic['Dimension']
variable = dic['Variable']

model_lst = ['Steady','PP','LM','O']
omega_lst = np.arange(0.05,0.301,0.05)
time_arr = np.arange(0,0.1,0.01) #Needs to be imported from results
timearr = time_arr = np.linspace(0,1,20) #Needs to be imported from results
var = ['alpha','phi','a','ap','f_tan','f_nor', 'CT', 'a_global']
labels = [r'$\alpha$ [deg]','$\phi$ [deg]', 'a [-]','$a^,[-]$', '$C_t$ [-]', '$C_n$ [-]', '$C_T$ [-]', 'a [-]']

for c in range(len(cases)):

    case_index=cases[c][0:3]

    if case_index == 'A.1':
        file = open("CT_step.pkl", "rb")
        data = pickle.load(file)
        magnitude_index = int(cases[0][4])-1
        mode = 'oscillations'


    elif case_index == 'A.2':
        file = open("CT_sin.pkl", "rb")
        data = pickle.load(file)
        freq_index = np.where(omega_lst == float(cases[0][6:]))[0][0]
        magnitude_index = int(cases[0][4])-1
        mode = 'oscillations'

    else:
        raise Exception('Incorrect case selected')

    for m in range(len(models)):

        model = models[m]
        model_index = np.where(np.asarray(model_lst) == model)[0][0]

        if case_index == 'A.1':
            #Ask Guillem to loop over multiple models too
            #result = data['results'][model_index][magnitude_index]
            result = data['results'][magnitude_index]

        elif case_index == 'A.2':
            result = data['results'][model_index][freq_index][magnitude_index]

        else:
            raise Exception('Third entry of case name probably wrong')

        if dimension[0]=='2D':
            plt.plot(time_arr, getattr(result, variable[0])) #Check non-dimensionalisation of time and output
            plt.xlabel('t [s]')
            lab = labels[np.where(np.asarray(var) == variable[0])[0][0]]
            plt.ylabel(lab)


        elif dimension[0] == '3D':
            Z = getattr(result, variable[0])
            Z = Z[:,0,:]

            fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
            mu = result.mu[:,0,0]
            r, t = np.meshgrid(time_arr, mu) #Check non-dimensionalisation of time
            ax.plot_surface(r, t, Z, cmap='viridis')

            ax.set_ylabel('$\mu$ [-]')
            ax.set_xlabel('t [s]')
            lab = labels[np.where(np.asarray(var) == variable[0])[0][0]]
            ax.set_zlabel(lab)

        else:
            raise Exception('Invalid dimension input')









