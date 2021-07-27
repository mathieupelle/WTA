
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.fft import fft, fftfreq
import pickle


x = 6
plt.rc('figure', figsize=[46.82 * .5**(.5 * x), 33.11 * .5**(.5 * x)])
plt.rc('font', family='serif')

save = True

def Plotting_following_blade(res,mod_labels,save_name,rad_pos = [0.4,0.6,0.8],lims=False):
    var = ['alpha','phi','a','ap','f_tan','f_nor','cl']
    labels = [r'$\alpha$ [deg]','$\phi$ [deg]', 'a [-]','$a^,[-]$', '$C_t$ [-]', '$C_n$ [-]', '$C_{l}$ [-]','$U_{rel}$ [m/s]']

    #cols = ['tab:blue','tab:orange','tab:green']
    cols = ['b','r','g']

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
            if par=='f_tan' or par=='f_nor':
                Z = Z/(0.5*1.225*10**2*50)
            x = np.zeros((len(rad_pos),len(res[j].time)))
            az_idx = 0
            for k in range(x.shape[1]): #We want to check how the blade position changes with time
                x[:,k] = Z[idx_mu.astype(int),az_idx,k]
                if az_idx+1 < Z.shape[1]:
                    az_idx = az_idx+1
                else:
                    az_idx = 0

            if lims:
                time = res[j].time/5-lims[0]

            else:
                time = res[j].time/5

            l0, = plt.plot(time,x[0,:],color=cols[j])
            line_handles.append(l0)
            try:
                l1, = plt.plot(time,x[1,:],'--',color=cols[j])
                try:

                    l2, = plt.plot(time,x[2,:],'-.',color=cols[j])
                    sec_leg_labels = list(np.round_(res[j].mu[idx_mu.astype(int),0,0],2))
                    sec_legend = plt.legend(handles=[l0,l1,l2], labels= sec_leg_labels)
                    # Add the legend manually to the current Axes.
                except:
                    sec_leg_labels = list(np.round_(res[j].mu[idx_mu.astype(int),0,0],2))
                    sec_legend = plt.legend(handles=[l0,l1], labels= sec_leg_labels,loc=2)
                    # Add the legend manually to the current Axes.
            except:
                pass


                # Create another legend for the second line.
                #plt.legend(handles=[line2], loc='lower right')
        plt.xlabel('Time [s]')
        plt.xlabel('Time $t U / R$ [-]')
        plt.ylabel(labels[i])
        plt.grid()
        plt.legend(handles = line_handles,labels = mod_labels,loc=1)
        try:
            plt.gca().add_artist(sec_legend)
        except:
            pass
        if lims:
            plt.xlim(list(np.array(lims)-lims[0]))
        plt.show()
        if save==True:
            plt.savefig('figures/'+save_name+'_'+par+'.pdf')




def Plotting_integral_quantities(res,mod_labels,save_name,lims=False):
    var = ['CT','CP']
    labels=[r'$C_T$',r'$C_P$']
    cols = ['b','r','g']

    for i,par in enumerate(var):
        plt.figure()
        for j in range(len(res)):
            Z =getattr(res[j],par)
            if lims:
                time = res[j].time/5-lims[0]
            else:
                time = res[j].time/5
            if j==5:
                plt.plot(time,Z,'--k', label=mod_labels[j],)
            else:
                plt.plot(time,Z,label=mod_labels[j],color=cols[j])

        plt.xlabel('Time $t U / R$ [-]')
        plt.ylabel(labels[i])
        plt.grid()
        plt.legend(loc='best')
        if lims:
            plt.xlim(list(np.array(lims)-lims[0]))
        plt.show()
        if save==True:
            plt.savefig('figures/'+save_name+'_'+par+'.pdf')


def Plotting_polars(res, mod_labels, rad_pos, omega,save_name):

    time = res[0].time
    T = 2*np.pi/omega
    idx1 = np.argmin(np.abs(time-1*T))-1
    #idx1 = 18
    idx2 = np.argmin(np.abs(time-2*T))+2
    #idx2 = len(time)
    cols = ['b','r','g']

    idx_mu = np.zeros(len(rad_pos))
    for i,val in enumerate(rad_pos):
        idx_mu[i] = (np.abs(res[0].mu[:,0,0] - val)).argmin()

    plt.figure()
    for j in range(len(res)):

        if j>0:
            Z = res[j].UA_class.Cn
        else:
            Z = res[j].cl

        x = np.zeros((len(rad_pos),len(res[j].time)))
        az_idx = 0
        for k in range(x.shape[1]):
            x[:,k] = Z[idx_mu.astype(int),az_idx,k]
            if az_idx+1 < Z.shape[1]:
                az_idx = az_idx+1
            else:
                az_idx = 0

        a = np.zeros((len(rad_pos),len(res[j].time)))
        az_idx = 0
        for k in range(a.shape[1]):
            a[:,k] = res[j].alpha[idx_mu.astype(int),az_idx,k]
            if az_idx+1 < res[j].alpha.shape[1]:
                az_idx = az_idx+1
            else:
                az_idx = 0

        for l in range(len(rad_pos)):
            plt.plot(a[l,idx1:idx2],x[l,idx1:idx2],'.-', label=mod_labels[j], color=cols[j])

    data_airfoil = pd.read_excel('polar DU95W180 (3).xlsx',header = 3,names=['alpha','Cl','Cd','Cm'])
    plt.plot(data_airfoil.alpha[21:], data_airfoil.Cl[21:], '--k', label='Experimental')
    #plt.xlim([a[:,idx1:idx2].min()-2,a[:,idx1:idx2].max()+2])
    plt.legend()
    plt.grid()
    plt.xlabel(r'$\alpha$ [deg]')
    plt.ylabel('$C_l$ [-]')
    if save==True:
        plt.savefig('figures/'+save_name+str(j)+'.pdf')

def Plotting_FFT(res, rad_pos, save_name):
    var = ['alpha','phi','a','ap','f_tan','f_nor','cl']
    cols = ['b','r','g']

    #Get the index of the wanted radial positions
    idx_mu = np.zeros(len(rad_pos))
    for i,val in enumerate(rad_pos):
        idx_mu[i] = (np.abs(res[0].mu[:,0,0] - val)).argmin()

    # Depends on the number of inputs that we pass
    for i,par in enumerate(var):
        plt.figure()
        for j in range(len(res)): #We want to compare all the results
            Z = getattr(res[j], par)
            if par=='f_tan' or par=='f_nor':
                Z = Z/(0.5*1.225*10**2*50)
            x = np.zeros((len(rad_pos),len(res[j].time)))
            az_idx = 0
            for k in range(x.shape[1]): #We want to check how the blade position changes with time
                x[:,k] = Z[idx_mu.astype(int),az_idx,k]
                if az_idx+1 < Z.shape[1]:
                    az_idx = az_idx+1
                else:
                    az_idx = 0

            for l in range(len(rad_pos)):
                yf = fft(x[l,:]-np.mean(x[l,:]))
                N = len(res[j].time)
                dt = res[j].time[1]
                xf = fftfreq(N, dt)*2*np.pi
                y_plot = []
                x_plot = []
                for q in range(len(yf)):
                    if 6>=(xf[q])>=0:
                        x_plot.append(xf[q])
                        y_plot.append(abs(yf[q])*2)

                plt.plot(x_plot, y_plot, color=cols[l], marker='.', markersize=10,label='$\mu$='+str(np.round(res[j].mu[idx_mu[l].astype(int),0,0],2)))


        plt.grid()
        plt.ylim(0)
        plt.xlim(0)
        plt.xlabel(r'$\omega$ [rad/s]')
        plt.ylabel(' Magnitude [s]')
        plt.legend()
        if save:
            plt.savefig('figures/'+save_name+'_'+par+'.pdf')



def Plotting_polars_radial(res, mod_labels, rad_pos, omega,save_name):

    time = res[0].time
    T = 2*np.pi/omega
    idx1 = np.argmin(np.abs(time-1*T))
    #idx1 = 18
    idx2 = np.argmin(np.abs(time-2*T))+1
    #idx2 = len(time)
    cols = ['b','g','r']

    plt.figure()
    idx_mu = np.zeros(len(rad_pos))
    for i,val in enumerate(rad_pos):
        idx_mu[i] = (np.abs(res[0].mu[:,0,0] - val)).argmin()


    for j in range(len(res)):

        if j>0:
            Z = res[j].UA_class.Cn
        else:
            Z = res[j].cl

        x = np.zeros((len(rad_pos),len(res[j].time)))
        az_idx = 0
        for k in range(x.shape[1]):
            x[:,k] = Z[idx_mu.astype(int),az_idx,k]
            if az_idx+1 < Z.shape[1]:
                az_idx = az_idx+1
            else:
                az_idx = 0

        a = np.zeros((len(rad_pos),len(res[j].time)))
        az_idx = 0
        for k in range(a.shape[1]):
            a[:,k] = res[j].alpha[idx_mu.astype(int),az_idx,k]
            if az_idx+1 < res[j].alpha.shape[1]:
                az_idx = az_idx+1
            else:
                az_idx = 0

        for l in range(len(rad_pos)):
            plt.plot(a[l,idx1:idx2],x[l,idx1:idx2],'.-', label=mod_labels[l], color=cols[l])

    data_airfoil = pd.read_excel('polar DU95W180 (3).xlsx',header = 3,names=['alpha','Cl','Cd','Cm'])
    plt.plot(data_airfoil.alpha[21:], data_airfoil.Cl[21:], '--k', label='Experimental')

    plt.legend()
    plt.grid()
    plt.xlabel(r'$\alpha$ [deg]')
    plt.ylabel('$C_l$ [-]')
    if save==True:
        plt.savefig('figures/'+save_name+'.pdf')


#%% Case 1 plotting

file = open("Dyn1.pkl", "rb")
data = pickle.load(file)
omega = 1.6
T = 2*np.pi/omega

Plotting_integral_quantities(data['results'],['Steady','DI', 'DS'],'org_noLims')

Plotting_following_blade(data['results'],['Steady', 'DS (noLEsep)','DS'],'org_Lims',rad_pos=[0.4, 0.8], lims=[T/5,2*T/5])

Plotting_polars(data['results'][0:3], ['Steady', 'DS (noLEsep)','DS'], rad_pos=[0.4], omega=omega, save_name='polar_comp')

Plotting_polars_radial(data['results'][2:3], ['$\mu$=0.4', '$\mu$=0.8'],rad_pos=[0.4,0.8], omega=omega, save_name='Dyn1_polars')



#%% Case 2 plotting

file = open("Dyn2.pkl", "rb")
data = pickle.load(file)
omega = 0.06
T = 2*np.pi/omega

Plotting_following_blade(data['results'][1:3],['DI','DS'],'Dyn2_rad0.2',rad_pos=[0.2], lims=[3*T/5,4*T/5])
Plotting_integral_quantities(data['results'],['Steady','DI', 'DS'],'Dyn2',lims=[3*T/5,4*T/5])
Plotting_polars(data['results'][1:3], ['DI', 'DS'],rad_pos=[0.], omega=0.06+1.7, save_name='polar')
Plotting_FFT(data['results'][2:3], [0.2,0.9], save_name='FFT')

