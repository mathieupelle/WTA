
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt



x = 6
plt.rc('figure', figsize=[46.82 * .5**(.5 * x), 33.11 * .5**(.5 * x)])
plt.rc('font', family='serif')


def Plotting_following_blade(res,mod_labels,rad_pos = [0.4,0.6,0.8],lims=False):
    var = ['alpha','phi','a','ap','f_tan','f_nor','cl']
    labels = [r'$\alpha$ [deg]','$\phi$ [deg]', 'a [-]','$a^,[-]$', '$C_t$ [-]', '$C_n$ [-]', '$C_{l}$ [-]']
    
    cols = ['tab:blue','tab:orange','tab:green']
        
    #Get the index of the wanted radial positions
    idx_mu = np.zeros(len(rad_pos))
    for i,val in enumerate(rad_pos):
        idx_mu[i] = (np.abs(res[0].mu[:,0,0] - val)).argmin()
    
    # Depends on the number of inputs that we pass
    for i,par in enumerate(var):
        plt.figure()
        line_handles = []
        for j in range(len(res)): #We want to compare all the results
            Z = getattr(res[j], par)
            x = np.zeros((len(rad_pos),len(res[j].time)))
            az_idx = 0
            for k in range(x.shape[1]): #We want to check how the blade position changes with time
                x[:,k] = Z[idx_mu.astype(int),az_idx,k]
                if az_idx+1 < Z.shape[1]:
                    az_idx = az_idx+1
                else:
                    az_idx = 0
            
            l0, = plt.plot(res[j].time,x[0,:],color=cols[j])
            line_handles.append(l0)
            try:
                l1, = plt.plot(res[j].time,x[1,:],'--',color=cols[j])
                try:
                    l2, = plt.plot(res[j].time,x[2,:],'-.',color=cols[j])   
                    sec_leg_labels = list(np.round_(res[j].mu[idx_mu.astype(int),0,0],2))
                    sec_legend = plt.legend(handles=[l0,l1,l2], labels= sec_leg_labels)
                    # Add the legend manually to the current Axes.
                except:
                    sec_leg_labels = list(np.round_(res[j].mu[idx_mu.astype(int),0,0],2))
                    sec_legend = plt.legend(handles=[l0,l1], labels= sec_leg_labels,loc=1)
                    # Add the legend manually to the current Axes.
            except:
                pass
            
            
                # Create another legend for the second line.
                #plt.legend(handles=[line2], loc='lower right')          
        plt.xlabel('Time [s]')
        plt.ylabel(labels[i])
        plt.grid()
        plt.legend(handles = line_handles,labels = mod_labels,loc=4)
        try:
            plt.gca().add_artist(sec_legend)
        except:
            pass
        if lims:
            plt.xlim(lims)
        plt.show()

        

def Plotting_integral_quantities(res,mod_labels,lims=False):
    var = ['CT','CP']
    labels=[r'$C_T$',r'$C_P$']
    
    for i,par in enumerate(var):
        plt.figure()
        for j in range(len(res)):
            Z =getattr(res[j],par)
            plt.plot(res[j].time,Z,label=mod_labels[j])
            
            
        plt.xlabel('Time [s]')
        plt.ylabel(labels[i])
        plt.grid()
        plt.legend(loc='best')
        if lims:
            plt.xlim(lims)
        plt.show()