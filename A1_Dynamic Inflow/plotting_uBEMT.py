# -*- coding: utf-8 -*-
"""
Created on Fri Jul  2 15:17:57 2021

@author: Mathieu Pellé
"""

import numpy as np
import matplotlib.pyplot as plt
import pickle

x = 6
plt.rc('figure', figsize=[46.82 * .5**(.5 * x), 33.11 * .5**(.5 * x)])
plt.rc('font', family='serif')

def plotting(plot_input,name, save=False):

    dic = plot_input

    cases = dic['Cases']
    models =  dic['Models']
    dimension = dic['Dimension']
    variable = dic['Variable']

    model_lst = ['Steady','PP','LM','O']
    omega_lst = np.array([0.05,0.15,0.3])
    var = ['alpha','phi','a','ap','f_tan','f_nor', 'CT', 'a_global']
    labels = [r'$\alpha$ [deg]','$\phi$ [deg]', 'a [-]','$a^,[-]$', '$C_t$ [-]', '$C_n$ [-]', '$C_T$ [-]', 'a [-]']
    models_legend = ['Steady', 'Pitt-Peters', 'Larsen-Madsen', 'Øye']
    step_legend = {'A.1.1': 'CT 0.5-->0.9','A.1.2': 'CT 0.9-->0.5', 'A.1.3': 'CT 0.2-->1.1', 'A.1.4': 'CT 1.1-->0.4',
                      'B.1.1': '$U_{\infty}/U_0$ 1-->1.5', 'B.1.2': '$U_{\infty}/U_0$ 1-->0.8',
                      'B.1.3': '$U_{\infty}/U_0$ 1-->1.2', 'B.1.4': '$U_{\infty}/U_0$ 1-->0.9'}
    oscillation_legend = {'A.2.1': ['0.5', '0.5'], 'A.2.2': ['0.9', '0.3'], 'A.2.3': ['0.2', '0.7'],
                          'B.2.1': ['1', '0.5'], 'B.2.2': ['0.8', '0.3'], 'B.2.3': ['1.2', '0.5']}

    if dimension[0]!='3D':
        plt.figure()
    else:
        pass
    for c in range(len(cases)):

        case_index=cases[c][0:3]

        if case_index == 'A.1':
            file = open("CT_step.pkl", "rb")
            data = pickle.load(file)
            magnitude_index = int(cases[c][4])-1
            mode = 'step'
            shift = 300
            print('CT Step')


        elif case_index == 'A.2':
            file = open("CT_sin.pkl", "rb")
            data = pickle.load(file)
            freq_index = np.where(omega_lst == float(cases[c][6:]))[0][0]
            magnitude_index = int(cases[c][4])-1
            mode = 'oscillations'
            print('CT Oscillations')

        elif case_index == 'B.1':
            file = open("U_step.pkl", "rb")
            data = pickle.load(file)
            magnitude_index = int(cases[c][4])-1
            #scaling = 1.1487626841233327**2
            scaling =  1
            mode = 'step'
            shift = len(data['time'])
            print('U step')


        elif case_index == 'B.2':
            file = open("U_sin.pkl", "rb")
            data = pickle.load(file)
            freq_index = np.where(omega_lst == float(cases[c][6:]))[0][0]
            magnitude_index = int(cases[c][4])-1
            scaling =  1
            mode = 'oscillations'
            print('U oscillations')


        else:
            raise Exception('Incorrect case selected')
        time_arr = data['time']*10/50
        if mode == 'step':
            time_arr = time_arr[:shift]
        for m in range(len(models)):
            model = models[m]
            model_index = np.where(np.asarray(model_lst) == model)[0][0]

            if case_index == 'A.1' or case_index=='B.1':
                result = data['results'][model_index][magnitude_index]

            elif case_index == 'A.2' or case_index=='B.2':
                result = data['results'][model_index][freq_index][magnitude_index]


            else:
                raise Exception('Third entry of case name probably wrong')

            if dimension[0]=='2D':
                print('2D plot')

                if variable[0]=='a_global':
                    Z = getattr(result, 'a')
                    Z = np.mean(Z,axis=0)[0,:]
                elif variable[0]=='CT':
                    Z = getattr(result, variable[0])
                    if case_index == 'B.1':
                        Z[1:] = Z[1:]*scaling
                else:
                    raise Exception('Not able to plot this variable in 2D time plot (yet?)')

                if len(cases)>1:
                    if mode=='step':
                        leg = step_legend[cases[c]]
                    elif mode=='oscillations':
                        lst = oscillation_legend[cases[c][0:5]]
                        if case_index=='A.2':
                            leg = 'CT='+str(lst[0])+'+'+str(lst[1])+'sin('+str(omega_lst[freq_index])+'t)'
                        else:
                            leg = '$U_{\infty}$='+str(lst[0])+'+'+str(lst[1])+'sin('+str(omega_lst[freq_index])+'t)'
                    else:
                        pass

                elif len(models)>1:
                    leg = models_legend[model_index]
                else:
                    pass
                if mode == 'step':
                    Z = Z[:shift]

                plt.plot(time_arr, Z, label=leg)
                plt.xlabel('$U_{\infty}t/R$ [-]')
                lab = labels[np.where(np.asarray(var) == variable[0])[0][0]]
                plt.ylabel(lab)


            elif dimension[0] == '3D':
                print('3D plot')
                Z = getattr(result, variable[0])
                Z = Z[:,0,:]
                if variable[0]=='f_tan' or variable[0]=='f_nor':
                    Z=Z/(0.5*1.225*10**2*50)

                _, ax = plt.subplots(subplot_kw={"projection": "3d"})
                mu = result.mu[:,0,0]
                if mode == 'step':
                    mu = mu[:shift]
                r, t = np.meshgrid(time_arr, mu)
                if mode == 'step':
                    Z = Z[:,:shift]
                #print(np.shape(r), np.shape(t), np.shape(Z))
                ax.plot_surface(r, t, Z, cmap='viridis')
                ax.set_ylim(ax.get_ylim()[::-1])

                ax.set_ylabel('$\mu$ [-]')
                ax.set_xlabel('$U_{\infty}t/R$ [-]')
                lab = labels[np.where(np.asarray(var) == variable[0])[0][0]]

                #ax.set_xlim(0,3)
                ax.view_init(elev=20., azim=-130)
                ax.zaxis.set_rotate_label(False)
                ax.set_zlabel(lab, rotation=0, labelpad=9)

            elif dimension[0]=='r_pos':
                print('2D radial position plot')

                if variable[0]=='a':
                    Z = getattr(result, 'a')
                    mu = getattr(result, 'mu')
                    idx = np.argmin(abs(mu[:,0,0]-dimension[1]))
                    Z = Z[idx,0,:]
                else:
                    raise Exception('Not able to plot this variable in 2D radial position plot (yet?)')

                if mode == 'step':
                    Z = Z[:shift]

                plt.plot(time_arr, Z, label=models_legend[model_index])
                plt.xlabel('$U_{\infty}t/R$ [-]')
                lab = labels[np.where(np.asarray(var) == variable[0])[0][0]]
                plt.ylabel(lab)

            else:
                raise Exception('Invalid dimension input')


    plt.grid()
    if dimension[0]!='3D':
        plt.legend()

    if dimension[0]=='2D':
        plt.xlim([0,max(time_arr)])
    if dimension[0]=='r_pos':
        plt.ylim([0,0.5])
        plt.xlim([0,max(time_arr)])
    if save:
        plt.savefig('figures/'+str(name)+'.pdf')


#%% Plotting cases definition

#Variables that can be plotted in 3D Z vs t vs r -->'alpha','phi','a','ap','f_tan','f_nor',
#Variables that can be plotted in 2D Z vs t -->'CT', 'a_global'

#Steps in CT
# A.1.1 0.5 --> 0.9
# A.1.2 0.9 --> 0.5
# A.1.3 0.2 --> 1.1
# A.1.4 1.1 --> 0.4

#Oscillations in CT
# A.x.x_omega --> omega = desired frequency
# A.2.1 --> CT0=0.5 dCT=0.5
# A.2.2 --> CT0=0.9 dCT=0.3
# A.2.3 --> CT0=0.2 dCT=0.7

#Steps in CT
# B.1.1 1 --> 1.5
# B.1.2 1 --> 0.8
# B.1.3 1 --> 1.2
# B.1.4 1 --> 0.9

#Oscillations in CT
# B.x.x_omega --> omega = desired frequency
# B.2.1 --> U1/U0=1 dU=0.5
# B.2.2 --> U1/U0=0.8 dU=0.3
# B.2.3 --> U1/U0=1.2 dU=0.5

saving=False

#%% Model comparison

plot_input = {'Cases':['A.1.1'], 'Models':['Steady','PP','LM','O'], 'Dimension':['2D'], 'Variable':['CT']}
plotting(plot_input,'CT_3models_CTstep1',save=saving)

plot_input = {'Cases':['A.1.1'], 'Models':['Steady','PP','LM','O'], 'Dimension':['2D'], 'Variable':['a_global']}
plotting(plot_input,'a_global_3models_CTstep1',save=saving)


plot_input = {'Cases':['A.1.1'], 'Models':['Steady','PP','LM','O'], 'Dimension':['r_pos',0.25], 'Variable':['a']}
plotting(plot_input,'a_global_0.25_3models_CTstep1',save=saving)

plot_input = {'Cases':['A.1.1'], 'Models':['Steady','PP','LM','O'], 'Dimension':['r_pos',0.5], 'Variable':['a']}
plotting(plot_input,'a_global_0.5_3models_CTstep1',save=saving)

plot_input = {'Cases':['A.1.1'], 'Models':['Steady','PP','LM','O'], 'Dimension':['r_pos',0.75], 'Variable':['a']}
plotting(plot_input,'a_global_0.75_3models_CTstep1',save=saving)

#%% Step comparison

plot_input = {'Cases':['A.1.1', 'A.1.2', 'A.1.3', 'A.1.4'], 'Models':['PP'], 'Dimension':['2D'], 'Variable':['CT']}
plotting(plot_input,'CT_PP_4CTsteps',save=saving)

plot_input = {'Cases':['A.1.1', 'A.1.2', 'A.1.3', 'A.1.4'], 'Models':['PP'], 'Dimension':['2D'], 'Variable':['a_global']}
plotting(plot_input,'a_global_PP_4CTsteps',save=saving)

v = ['alpha','phi','a','ap','f_tan','f_nor']
for i in range(len(v)):
    plot_input = {'Cases':['A.1.1'], 'Models':['PP'], 'Dimension':['3D'], 'Variable':[v[i]]}
    plotting(plot_input,v[i]+'_PP_3D_CTstep1',save=saving)

#%% Oscillation comparison

plot_input = {'Cases':['A.2.2_0.05', 'A.2.2_0.15', 'A.2.2_0.3'], 'Models':['PP'], 'Dimension':['2D'], 'Variable':['CT']}
plotting(plot_input,'CT_PP_3freq_amp1',save=saving)

plot_input = {'Cases':['A.2.2_0.05', 'A.2.2_0.15', 'A.2.2_0.3'], 'Models':['PP'], 'Dimension':['2D'], 'Variable':['a_global']}
plotting(plot_input,'a_global_PP_3freq_amp0.5',save=saving)

v = ['alpha','phi','a','ap','f_tan','f_nor']
for i in range(len(v)):
    plot_input = {'Cases':['A.2.1_0.3'], 'Models':['PP'], 'Dimension':['3D'], 'Variable':[v[i]]}
    plotting(plot_input,v[i]+'_PP_3D_freq0.3_amp0.5',save=saving)


#%% U_inf comparison

plot_input = {'Cases':['B.1.1'], 'Models':['Steady','PP','LM','O'], 'Dimension':['2D'], 'Variable':['CT']}
plotting(plot_input,'U_3models_Ustep1',save=saving)

plot_input = {'Cases':['B.1.1'], 'Models':['Steady','PP','LM','O'], 'Dimension':['2D'], 'Variable':['a_global']}
plotting(plot_input,'U_3models_Ustep1',save=saving)

plot_input = {'Cases':['B.1.1','B.1.2', 'B.1.4', 'B.1.4'], 'Models':['PP'], 'Dimension':['2D'], 'Variable':['CT']}
plotting(plot_input,'U_PP_4CTsteps',save=saving)

plot_input = {'Cases':['B.2.1_0.05', 'B.2.1_0.15', 'B.2.1_0.3'], 'Models':['PP'], 'Dimension':['2D'], 'Variable':['CT']}
plotting(plot_input,'U_PP_4CTsteps',save=saving)

plot_input = {'Cases':['B.2.1_0.05', 'B.2.1_0.15', 'B.2.1_0.3'], 'Models':['PP'], 'Dimension':['2D'], 'Variable':['a_global']}
plotting(plot_input,'U_PP_4CTsteps',save=saving)

