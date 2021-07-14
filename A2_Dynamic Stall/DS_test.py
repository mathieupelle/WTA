
import numpy as np
import matplotlib.pyplot as plt
from utilities_uBEMT import Rotor,BEMT
import pickle


Geometry = Rotor(N_radial_sections = 15)
Calc = BEMT(Geometry)

Calc.CpLambda(TSR_list = list(np.linspace(5,15,20)), theta_list = list(np.linspace(-7,5,20)))

#%%

U_0 = 10
TSR_0 = 8
omega = 0.3
T = 2*np.pi/omega
time_arr = np.arange(0,2*T,0.1)
Omega_0 = TSR_0*U_0/Geometry.radius

#Calculate the pitch angle corresponding to the optimal CT=8/9
pitch =  Calc.getPitchAngle_fromCT(CT = 8/9,TSR = TSR_0)
#U_arr = (1 + 0.5*np.cos(omega*time_arr))*U_0

phi = np.linspace(0,2*np.pi,9)
U_arr = np.zeros((8, len(time_arr)))
for i in range(len(phi)-1):
    azimuth = (phi[i]+phi[i+1])/2
    U_arr[i,:] = (1 + 0.5*np.cos(omega*time_arr)*np.cos(azimuth))*U_0


TSR_arr = Omega_0*Geometry.radius/(1 + 0.5*np.cos(omega*time_arr)*U_0)
#Build the conditions dict necessary for the unsteady BEMT
cond = {'wind_speed': U_arr,
        'pitch_angle': pitch*np.ones(len(time_arr)),
        'yaw_angle': np.zeros(len(time_arr)),
        'TSR': TSR_arr,
        'omega': Omega_0*np.ones(len(time_arr))}

Calc.Solver(time = time_arr, conditions = cond, DI_Model = 'PP', DS_Model='BL')

U_sin_results = Calc.Results
file = open("DS_oscU2.pkl","wb")
pickle.dump(U_sin_results,file)
file.close()
#%%

x = 6
plt.rc('figure', figsize=[46.82 * .5**(.5 * x), 33.11 * .5**(.5 * x)])
plt.rc('font', family='serif')


# Ncycles = np.floor(time_arr[-1]*omega/(2*np.pi))
# n_of_cycle = time_arr*omega/(2*np.pi)
# i1=np.argmin(np.abs(n_of_cycle-(Ncycles-1)))
# i2=np.argmin(np.abs(n_of_cycle-(Ncycles-.5)))
# i3=np.argmin(np.abs(n_of_cycle-(Ncycles)))+1

omega = 0.3
T = 2*np.pi/omega
i1 = np.argmin(np.abs(time_arr-1*T))
i3 = np.argmin(np.abs(time_arr-2*T))
# i1=0
# i3=len(time_arr)

result = U_sin_results
#file = open("DS_constantU.pkl", "rb")
#file = open("DS_constantU_TSR10.pkl", "rb")
#file = open("DS_oscU.pkl", "rb")
#file = open("DS_oscU_TSR10.pkl", "rb")
#file = open("DS_oscU_noLEsep.pkl", "rb")
#result = pickle.load(file)

plt.figure()
for i in range(len(time_arr)):
    if i1<=i<=i3:
        if np.max(result.ite[:,0,i])>500:
            plt.plot(np.linspace(0,1,14), result.ite[:,0,i], label = i)
        else:
            plt.plot(np.linspace(0,1,14), result.ite[:,0,i])
    else:
        pass
plt.grid()
plt.xlabel('r/R [-]')
plt.ylabel('Iterations [-]')
plt.legend()

var = ['alpha','phi','a','ap','f_tan','f_nor','cl']
labels = [r'$\alpha$ [deg]','$\phi$ [deg]', 'a [-]','$a^,[-]$', '$C_t$ [-]', '$C_n$ [-]', '$C_{l}$ [-]']

for i in range(len(var)):

    _, ax = plt.subplots(subplot_kw={"projection": "3d"})
    mu = result.mu[:,0,0]
    r, t = np.meshgrid(time_arr[i1:i3], mu)
    Z = getattr(result, var[i])
    Z = Z[:,3,i1:i3]

    ax.plot_surface(r, t, Z, cmap='viridis')
    #ax.set_ylim(ax.get_ylim()[::-1])

    ax.set_ylabel('$\mu$ [-]')
    ax.set_xlabel('$t$ [s]')
    ax.view_init(elev=20., azim=-130)
    ax.zaxis.set_rotate_label(False)
    ax.set_zlabel(labels[i], rotation=0, labelpad=9)

#%%
import pandas as pd

data_airfoil = pd.read_excel('polar DU95W180 (3).xlsx',header = 3,names=['alpha','Cl','Cd','Cm'])
plt.figure()
plt.plot(data_airfoil.alpha[21:], data_airfoil.Cl[21:], '--k')
alpha = getattr(result,'alpha')
unsteady_polar = getattr(result,'UA_class')
for i in range(14):
    plt.plot(alpha[i,0,i1:i3], unsteady_polar.Cn[i,15,i1:i3], label=round(mu[i],2))
plt.grid()
plt.xlabel(r'$\alpha$ [deg]')
plt.ylabel('$C_{l}$ [-]')
plt.legend()





