# -*- coding: utf-8 -*-
"""
Created on Thu Jul  8 16:03:22 2021

@author: Mathieu Pellé
"""

import numpy as np

class UnsteadyAerodynamics:
    def __init__(self, time, N_radial):
        """
          Class that computes unsteady airfoil performance.

          Parameters
          ----------
          time : Time vector
          N_radial: Number of radial elements

        """

        self.time = time
        #Empty arrays to store parameters
        self.X, self.Y, self.dalpha_qs_dt, self.D, self.delta_dalpha_qs_dt, self.Dp, self.Cn_P, self.f, self.Dbl, self.tau_v, self.Cvortex, self.Cn_vortex, self.Cn = np.zeros((13, N_radial, len(self.time)+1))

    def Duhamel(self, X, Y, ds, dalpha, order=2, A1=0.3, b1=0.14, A2=0.7, b2=0.53):
        """
        Duhamel's superposition integral with Wagner's indicial function response to step changes in AoA.

        Parameters
        ----------
        X : Previous lag state X.
        Y : Previous lag state Y.
        ds : Semichord step.
        dalpha : Quasi-steady alpha step.
        order : Order of equation. The default is 2.
        A1 : Coefficient for X state. The default is 0.3.
        b1 : Exponent for X state. The default is 0.14.
        A2 : Coefficient for Y state. The default is 0.7.
        b2 : Exponent for Y state The default is 0.53.

        Returns
        -------
        X : New lag state X.
        Y : New lag state Y.

        """

        # A1 = 0.165
        # A2 = 0.335
        # b1 = 0.335
        # b2 = 0.3

        if order==1:
            X = X*np.exp(-b1*ds)+A1*dalpha
            Y = Y*np.exp(-b2*ds)+A2*dalpha

        elif order==2:

            X = X*np.exp(-b1*ds)+A1*dalpha*np.exp(-b1*ds/2)
            Y = Y*np.exp(-b2*ds)+A2*dalpha*np.exp(-b2*ds/2)
        else:

            X = X*np.exp(-b1*ds)+A1*dalpha*((1+4*np.exp(-b1*ds/2)+np.exp(-b1*ds))/6)
            Y = Y*np.exp(-b2*ds)+A2*dalpha*((1+4*np.exp(-b2*ds/2)+np.exp(-b2*ds))/6)

        return X, Y


    #NC terms
    def Deficiency_NC(self, D, delta_dalpha_qs_dt, dt, c, Kalpha=0.75):
        """
        Computes value of deficiency function used fo non-circulatory loads.

        Parameters
        ----------
        D : Previous deficiency term.
        delta_dalpha_qs_dt : Difference between quasi steady AoA gradients at current and previous time steps.
        dt : Time step
        c : Chord length
        Kalpha : Factor. The default is 0.75.

        Returns
        -------
        D : New deficiency term.
        Kalpha : Factor.

        """
        a = 343 #speed of sound
        Kalpha = 0.75
        D = D*np.exp(-a*dt/Kalpha/c)+delta_dalpha_qs_dt*np.exp(-a*dt/Kalpha/c/2)

        return D, Kalpha



    def f_sep(self, alpha, parameters = False):
        """
        Computes pseudo location of seperation poinit from angle of attack.

        Parameters
        ----------
        alpha : Angle of attack in radians.
        parameters : Dictionnary of tunning parameters. The default is False.

        Returns
        -------
        f : Pseudo-location of seperation.

        """
        # a defines angle of each parts of the curve
        # n affects gradient of curve
        # c affects how much effect the correction has
        if parameters:
            a1, a2, a3 = parameters['a1'], parameters['a2'], parameters['a3']
            n2, n3 = parameters['n2'], parameters['n3']
            c2, c3, c4 = parameters['c2'], parameters['c3'], parameters['c4']

        else:
            a1, a2, a3 = 7, 15, 21
            n2, n3 = 1, 0.3
            c2, c3, c4 = 0.8, 1, 0.0

        alpha = np.rad2deg(alpha)
        if alpha<=a1:
            f = 1
        elif alpha>a1 and alpha<=a2:
            f = 1-c2*((alpha-a1)/(a2-a1))**n2
        elif alpha>a2 and alpha<a3:
            f = 0.2*(1-c3*((alpha-a2)/(a3-a2))**n3)
        else:
            f = c4

        return f



    def Pressure_lag(self, Dp, ds, delta_Cn, Tp=1.7):
        """
        Computes deficiency for leading-edge pressure lag.

        Parameters
        ----------
        Dp : Previous deficiency term.
        ds : Semichord step.
        delta_Cn : Difference in normal coefficient at current and previous times.
        Tp : Factor dependent on Mach number and airfoil shape. The default is 1.7.

        Returns
        -------
        Dp : New deficiency term.

        """

        Dp = Dp*np.exp(-ds/Tp)+delta_Cn*np.exp(-ds/2/Tp)

        return Dp


    def BL_lag(self, Dbl, ds, delta_f, Tf = 3.0):
        """
        Computes deficiency for boundary-layer lag.

        Parameters
        ----------
        Dbl : Previous deficiency term.
        ds : Semichord step.
        delta_f : Difference in pseudo location of seperation point at current and previous times.
        Tf : Factor dependent on airfoil shape. The default is 3.0.

        Returns
        -------
        Dbl : New deficiency term.

        """

        Dbl = Dbl*np.exp(-ds/Tf)+delta_f*np.exp(-ds/2/Tf)

        return Dbl

    def Vortex_time(self, tau_v, Cn, ds, delta_alpha_qs, U_conv=0.45, Cn1=1.0093):
        """
        Parameters
        ----------
        tau_v : Previous reduced time term.
        Cn : Current normal coefficient.
        ds : Semichord step
        delta_alpha_qs : Difference between quasi steady AoA gradients at current and previous time steps.
        U_conv : Wake vortex convection speed. The default is 0.45.
        Cn1 : Critical normal coefficient. The default is 1.0093.

        Returns
        -------
        tau_v : New reduced time term.

        """
        if Cn>Cn1:
            tau_v = tau_v+U_conv*ds
        else:
            if delta_alpha_qs<0 and tau_v>0:
                tau_v = tau_v+U_conv*ds
            else:
                tau_v = 0

        return tau_v


    def Vortex_shedding(self, Cn, ds, delta_Cn, tau_v, TVL=11, TV=6):
        """
        Computes added loading from leading-edge vortex.

        Parameters
        ----------
        Cn : Current normal coefficient.
        ds : Semichord step
        delta_Cn : Difference in forcing term between current and previous times.
        tau_v : Current reduced time.
        TVL : Some kind of time limit? The default is 11.
        TV : Factor. The default is 6.

        Returns
        -------
        Cn : Normal coefficient from vortex.

        """

        if tau_v>0.001 and tau_v<TVL:
            Cn = Cn*np.exp(-ds/TV)+delta_Cn*np.exp(-ds/2/TV)

        else:
            Cn = Cn*np.exp(-ds/TV)

        return Cn



    def Beddoes_Leishman(self, airfoil, i, t, alpha, h, Uinf, LE_sep=True):
        """
        Applies the Beddoes_Leishman dynamic stall model to specified airfoil and for specified conditions,

        Parameters
        ----------
        t : Time index
        alpha : Angle of attack array.
        h : Heave displacement array.
        LE_sep : (De)Activate leading-edge seperation module. The default is True.

        Returns
        -------
        None.

        """
        #Extract airfoil parameters
        dCn_dalpha = airfoil['dCn_dalpha']
        c= airfoil['chord']
        alpha0 = np.deg2rad(airfoil['alpha0'])

        #Airfoil tuning
        parameters = {'a1':9.9, 'a2':17, 'a3':23.05,
                      'n2':0.45, 'n3':0.6,
                      'c2': 0.82, 'c3': 0.75, 'c4':0.01}

        if t==0:
            dt = self.time[t+1]-self.time[t] #Time step
            dalpha_dt = (alpha[t]-0)/dt #AoA derivative
            dh_dt = (h[t]-0)/dt #Heave derivative
        else:
            dt = self.time[t]-self.time[t-1] #Time step
            dalpha_dt = (alpha[t]-alpha[t-1])/dt #AoA derivative
            dh_dt = (h[t]-h[t-1])/dt #Heave derivative

        ds = 2*Uinf*dt/c #Time to semichords

        #Quasi-steayd AoA from airfoil motion
        alpha_qs1 = alpha[t-1] + dalpha_dt*c/2/Uinf-dh_dt/Uinf
        alpha_qs2 = alpha[t] + dalpha_dt*c/2/Uinf-dh_dt/Uinf
        self.dalpha_qs_dt[i, t] = (alpha_qs2-alpha_qs1)/dt

        #### Unsteady attached flow ####
        self.X[i, t], self.Y[i, t] = self.Duhamel(self.X[i, t-1], self.Y[i, t-1], ds, alpha_qs2-alpha_qs1) #Lag states
        delta_dalpha_qs_dt = self.dalpha_qs_dt[i, t]- self.dalpha_qs_dt[i, t-1]
        self.D[i, t], Kalpha = self.Deficiency_NC(self.D[i, t-1], delta_dalpha_qs_dt, dt, c) #Deficiency for NC loads

        alpha_eq = alpha_qs2 - self.X[i, t] - self.Y[i, t] #Equivalent AoA (lag from wake effects)
        Cn_C = dCn_dalpha*(alpha_eq-alpha0) #Circulatory loads
        Cn_NC = 4*Kalpha*c/Uinf*(self.dalpha_qs_dt[i, t]-self.D[i, t]) #Non-circulatory loads
        self.Cn_P[i, t] = Cn_NC + Cn_C #Total unsteady loads

        #### Non-linear TE seperation ####
        self.Dp[i, t] = self.Pressure_lag(self.Dp[i, t-1], ds, self.Cn_P[i, t]-self.Cn_P[i, t-1]) #Pressure lag deficiency

        alpha_f = (self.Cn_P[i, t]-self.Dp[i, t])/dCn_dalpha+alpha0
        self.f[i, t] = self.f_sep(alpha_f, parameters=parameters) #Seperation point

        self.Dbl[i, t] = self.BL_lag(self.Dbl[i, t-1], ds, self.f[i, t]-self.f[i, t-1]) #Boundary layer lag deficiency

        f_TE = self.f[i, t]-self.Dbl[i, t] #New seperation position

        Cn_f = dCn_dalpha*((1+np.sqrt(f_TE))/2)**2*(alpha_eq-alpha0)+Cn_NC #Total unsteady loads with TE seperation


        if LE_sep:

            #### LE seperation ####
            self.tau_v[i, t] = self.Vortex_time(self.tau_v[i, t-1], self.Cn_P[i, t-1]-self.Dp[i, t-1], ds, self.dalpha_qs_dt[i, t]-self.dalpha_qs_dt[i, t-1]) #??? #Reduced time

            self.Cvortex[i, t] = Cn_C*(1-((1+np.sqrt(f_TE))/2)**2) #Forcing term
            if t==0:
                self.Cn_vortex[i, t] = self.Cvortex[i, t]

            #### Vortex shedding ####
            self.Cn_vortex[i, t] = self.Vortex_shedding(self.Cn_vortex[i, t-1], ds, self.Cvortex[i, t]-self.Cvortex[i, t-1], self.tau_v[i, t-1]) #Vortex lift
            self.Cn[i, t] = Cn_f+self.Cn_vortex[i, t] #Unsteady loads with TE and LE seperations

        else:
            self.Cn[i, t] = Cn_f #Unsteady loads with TE seperation




#%%

import matplotlib.pyplot as plt
import pandas as pd

c = 1
airfoil = {'dCn_dalpha':2*np.pi, 'chord':c, 'alpha0':-2}
time = np.arange(0,500,0.1)
Uinf = 1

k = 0.1
omega = 2*k*c/Uinf
alpha = np.deg2rad(15)+np.deg2rad(10)*np.sin(omega*time)

k_h = 0.0
omega_h = 0
h = 0.3*np.sin(omega_h*time)


unsteady_polar = UnsteadyAerodynamics(time, 1)

for i in range(1):
    airfoil = {'dCn_dalpha':2*np.pi, 'chord':c, 'alpha0':-2}
    for t in range(len(time)):
        unsteady_polar.Beddoes_Leishman(airfoil, i, t, alpha, h, Uinf, LE_sep=True)

Ncycles = np.floor(time[-1]*omega/(2*np.pi))
n_of_cycle = time*omega/(2*np.pi)
i1=np.argmin(np.abs(n_of_cycle-(Ncycles-1)))
i2=np.argmin(np.abs(n_of_cycle-(Ncycles-.5)))
i3=np.argmin(np.abs(n_of_cycle-(Ncycles)))+1

data_airfoil = pd.read_excel('polar DU95W180 (3).xlsx',header = 3,names=['alpha','Cl','Cd','Cm'])
plt.figure()
plt.plot(data_airfoil.alpha[21:], data_airfoil.Cl[21:], '--k')
plt.plot(np.rad2deg(alpha)[i1:i3], unsteady_polar.Cn[0,i1:i3],'-b')
plt.grid()
plt.xlabel(r'$\alpha$ [deg]')
plt.ylabel('$C_{l}$ [-]')






