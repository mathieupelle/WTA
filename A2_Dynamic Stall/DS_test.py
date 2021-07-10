
import numpy as np
import matplotlib.pyplot as plt
from utilities_uBEMT import Rotor,BEMT


Geometry = Rotor(N_radial_sections = 15)
Calc = BEMT(Geometry)

Calc.CpLambda(TSR_list = list(np.linspace(5,15,20)), theta_list = list(np.linspace(-7,5,20)))

#%%

U_0 = 10
TSR_0 = 8
omega = 0.3
time_arr = np.arange(0,500,0.1)
Omega_0 = TSR_0*U_0/Geometry.radius

#Calculate the pitch angle corresponding to the optimal CT=8/9
pitch =  Calc.getPitchAngle_fromCT(CT = 8/9,TSR = TSR_0)
U_arr = (1 +0.5*np.sin(omega*time_arr))*U_0
TSR_arr = Omega_0*Geometry.radius/U_arr
#Build the conditions dict necessary for the unsteady BEMT
cond = {'wind_speed': U_arr,
        'pitch_angle': pitch*np.ones(len(time_arr)),
        'yaw_angle': np.zeros(len(time_arr)),
        'TSR': TSR_arr}

Calc.Solver(time = time_arr, conditions = cond, DI_Model = 'PP', DS_Model='BL')

U_sin_resuts = Calc.Results

#%%

x = 6
plt.rc('figure', figsize=[46.82 * .5**(.5 * x), 33.11 * .5**(.5 * x)])
plt.rc('font', family='serif')


# plt.figure()
# for i in range(len(time_arr)):
#     plt.plot(np.linspace(0,1,14), U_sin_resuts.ite[:,0,i])
# plt.grid()
# plt.xlabel('r/R [-]')
# plt.ylabel('Iterations [-]')

# plt.figure()
# for i in range(len(time_arr)):
#     plt.plot(np.linspace(0,1,14), U_sin_resuts.a[:,0,i])
# plt.grid()
# plt.xlabel('r/R [-]')
# plt.ylabel('a [-]')

# plt.figure()
# for i in range(len(time_arr)):
#     plt.plot(np.linspace(0,1,14), U_sin_resuts.f_nor[:,0,i])
# plt.grid()
# plt.xlabel('r/R [-]')
# plt.ylabel('Cn [-]')

Ncycles = np.floor(time_arr[-1]*omega/(2*np.pi))
n_of_cycle = time_arr*omega/(2*np.pi)
i1=np.argmin(np.abs(n_of_cycle-(Ncycles-1)))
i2=np.argmin(np.abs(n_of_cycle-(Ncycles-.5)))
i3=np.argmin(np.abs(n_of_cycle-(Ncycles)))+1

result = U_sin_resuts

var = ['alpha','phi','a','ap','f_tan','f_nor']
labels = [r'$\alpha$ [deg]','$\phi$ [deg]', 'a [-]','$a^,[-]$', '$C_t$ [-]', '$C_n$ [-]', '$C_T$ [-]', 'a [-]']

for i in range(len(var)):

    _, ax = plt.subplots(figsize=(8,5),subplot_kw={"projection": "3d"})
    mu = result.mu[:,0,0]
    r, t = np.meshgrid(time_arr[i1:i3], mu)
    Z = getattr(result, var[i])
    Z = Z[:,0,i1:i3]

    ax.plot_surface(r, t, Z, cmap='viridis')
    ax.set_ylim(ax.get_ylim()[::-1])

    ax.set_ylabel('$\mu$ [-]')
    ax.set_xlabel('$t$ [s]')
    ax.view_init(elev=20., azim=-130)
    ax.zaxis.set_rotate_label(False)
    ax.set_zlabel(labels[i], rotation=0, labelpad=9)
