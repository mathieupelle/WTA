

# -*- coding: utf-8 -*-
"""
Created on Sat Jun  5 11:19:26 2021

@author: MathieuPelle
"""

import numpy as np
import pandas as pd
import math as m
import warnings
import matplotlib.pyplot as plt


class rotor:
  def __init__(self, N_radial_sections = 50, Spacing_method = 'lin', pitch = -2):
    """
      Class that defines the geometry of a rotor.

      Parameters
      ----------
      Optimized_geometry : Optimal geometry determined using the class Optimizer.
          If input provided (none by default), the rotor will be created with that geometry
      N_radial_sections : Number of radial elements to be used
      Spacing method : 'Lin' -> linear spacing. 'cos' -> cosinusoidal spacing

   """

          #Geometric data of the blade
    self.radius = 50 #[m]
    self.n_blades = 3
    self.theta = pitch #Pitch angle [deg]
    self.N_radial = N_radial_sections #Number of sections

    #Create non-dimensional radius array depending on spacing method
    self.mu = np.linspace(0.2,1,self.N_radial)
    if Spacing_method == 'cos':
        angle=np.linspace(0,np.pi,self.N_radial)
        for i in range(len(self.mu)):
            self.mu[i]=0.4*(1-m.cos(angle[i]))+0.2
    self.beta = 14*(1-self.mu) #Twist angle in degrees
    self.chord = 3*(1-self.mu)+1 #Chord length in meters

    self.N_azimuth = 40 #Number of angular sections
    self.azimuth = np.linspace(0,2*np.pi,self.N_azimuth)

    #Polar data
    self.polars = pd.read_excel('polar DU95W180 (3).xlsx',header = 3,names=['alpha','Cl','Cd','Cm'])


    self.SetOperationalData(wind_speed=10, TSR=10) #Assign default values to operational conditions

  def SetOperationalData(self, wind_speed, TSR, rho=1.225):
    """
      Operational data associated to the rotor

      Parameters
      ----------
      wind_speed : Float [m/s]
      TSR : Float [-]
      rho : Float [kg/m3], optional
           The default is 1.225.
    """
    self.wind_speed = wind_speed
    self.TSR = TSR
    self.omega = wind_speed*TSR/self.radius
    self.rho = rho

class results: #Create the variables to store the results from BEMT
    def __init__(self, T, N_radial):
        self.a,self.ap,self.phi,self.alpha,self.cl,self.cd,self.f_nor,self.f_tan,self.f,self.f_tip,self.f_root,self.ite,self.chord,self.beta,self.mu,self.circulation,self.enthalpy_3,self.local_CT,self.local_CQ =  np.zeros((19,T,N_radial-1))

    def Integrate(self,Rotor_lst, T):
        #Calculate global CT
        #Generate differential radius array (to take into account cosinusoidal spacing)
        Rotor=Rotor_lst[0]
        d_r = np.zeros(len(Rotor.mu)-1)
        for i in range(len(Rotor.mu)-1):
            d_r[i] = (Rotor.mu[i+1]-Rotor.mu[i])*Rotor.radius

        self.CT, self.CP, self.CQ = np.zeros((3, T,))
        for t in range(T):
            Rotor = Rotor_lst[t]
            #print(Rotor.wind_speed)
            self.CT[t] = np.sum(self.f_nor[t,:].transpose()*Rotor.n_blades*d_r)/(0.5*Rotor.rho*Rotor.wind_speed**2*np.pi*Rotor.radius**2)

            #Global CP
            dTorque = self.f_tan[t,:].transpose()*d_r*self.mu[t,:]*Rotor.radius

            self.CP[t] = np.sum(dTorque*Rotor.n_blades*Rotor.omega)/(0.5*Rotor.rho*Rotor.wind_speed**3*np.pi*Rotor.radius**2)

            #Global CQ
            self.CQ[t] = np.sum(dTorque*Rotor.n_blades)/(0.5*Rotor.rho*Rotor.wind_speed**2*np.pi*Rotor.radius**3)


def AirfoilCoefficients(Rotor, alpha):
    cl = np.interp(alpha*180/np.pi,np.array(Rotor.polars['alpha']),np.array(Rotor.polars['Cl']))
    cd = np.interp(alpha*180/np.pi,np.array(Rotor.polars['alpha']),np.array(Rotor.polars['Cd']))

    return cl, cd


def RelativeVelocities(Rotor, a, ap, r):
    u_nor = Rotor.wind_speed*(1-a)
    u_tan = Rotor.omega*r*(1+ap)
    u_rel = m.sqrt(u_nor**2+u_tan**2)
    phi = m.atan(u_nor/u_tan)

    return u_tan, u_nor, u_rel, phi

def Forces(Rotor, chord, phi, u_rel, cl, cd):
    lift = 0.5*Rotor.rho*u_rel**2*chord*cl
    drag = 0.5*Rotor.rho*u_rel**2*chord*cd

    f_tan = lift*np.sin(phi) - drag*np.cos(phi)
    f_nor = lift*np.cos(phi) + drag*np.sin(phi)

    return lift, drag, f_tan, f_nor

def PrandtlTipCorrection(Rotor, mu, a_new):
    mu_root = Rotor.mu[0]
    #Tip correction
    exp = np.exp(-Rotor.n_blades/2 * ((1-mu)/mu) * np.sqrt(1+Rotor.TSR**2*mu**2/(1-a_new)**2))
    f_tip = 2/np.pi * np.arccos(exp)
    #Root correction
    exp = np.exp(-Rotor.n_blades/2 * ((mu-mu_root)/mu) * np.sqrt(1+Rotor.TSR**2*mu**2/(1-a_new)**2))
    f_root = 2/np.pi * np.arccos(exp)
    #Combined correction
    f = f_tip*f_root

    if f < 1e-4 or m.isnan(f):
        f = 1e-4

    return f,f_tip,f_root

def NewInductionFactor(CT):
    CT_1 = 1.816
    CT_2 = 2*np.sqrt(CT_1) - CT_1
    if CT < CT_2:
        a_new = 0.5 - np.sqrt(1-CT)/2
    else:
        a_new = 1 + (CT-CT_1)/(4*np.sqrt(CT_1)-4)
    return a_new

def NewCT(a, glauert=True):
    CT = 4*a*(1-a)
    if glauert:
        CT_1 = 1.816
        a_1 = 1-np.sqrt(CT_1)/2
        if a>a_1:
            CT = CT_1-4*(np.sqrt(CT_1)-1)*(1-a)
    return CT


def Pitt_Peters(CT, a, r, dt, Rotor):
    CT_qs = NewCT(a)
    dvdt = 3*m.pi*Rotor.wind_speed**2/16/r*(CT - CT_qs)
    #print(dvdt)
    v = -a*Rotor.wind_speed - dvdt*dt
    a_new = -v/Rotor.wind_speed
    return a_new

def plot_TSR(results, param_lst, times, save=False):
    var=['alpha','phi','a','ap','f_tan','f_nor','circulation','local_CQ','local_CT']
    labels=[r'$\alpha$ [deg]','$\phi$ [deg]', 'a [-]','$a^,[-]$', '$C_t$ [-]', '$C_n$ [-]','$\Gamma$ [-]','$C_q [-]$', '$C_T [-]$']
    for i in range(len(var)):
        plt.figure()
        plt.grid()
        plt.xlabel(r'Radius $\frac{r}{R}$ [-]')
        plt.ylabel(labels[i])
        for j in range(len(times)):
            dic=results
            param = param_lst[times[j]]
            if var[i]=='f_tan' or var[i]=='f_nor':
                Z=getattr(dic, str(var[i]))/(0.5*param.rho*param.wind_speed**2*param.radius)
            elif var[i]=='circulation':
                Z=getattr(dic, str(var[i]))/((np.pi*param.wind_speed**2/(param.n_blades*param.omega)))
            else:
                Z=getattr(dic, str(var[i]))
            plt.plot(dic.mu[times[j],:],Z[times[j],:],label='$t$=' +str(times[j]))

        plt.legend()
        if save==True:
            plt.savefig('figures/TSR_'+str(var[i])+'.pdf')
#%%

N_radial = 40

Rotor_lst = []
TSR = 10
wind_speed = 10
Res = {}

dt = 0.1
time = np.arange(0, 0.4, dt)

wind_vector = wind_speed*np.ones(len(time))
#wind_vector[1:] = wind_speed*1.5

pitch_lst = -2*np.ones(len(time))
pitch_lst[1:] = 0


Results = results(len(time), N_radial)
for t in range(len(time)):
    Rotor = rotor(N_radial_sections=N_radial, pitch=pitch_lst[t])
    Rotor.SetOperationalData(wind_vector[t],TSR)
    Rotor_lst.append(Rotor)


    for i in range(Rotor.N_radial-1):
        mu = (Rotor.mu[i]+Rotor.mu[i+1])/2
        r = Rotor.radius*mu
        chord = np.interp(mu,Rotor.mu,Rotor.chord)
        beta = np.interp(mu,Rotor.mu,Rotor.beta)
        [Results.mu[t, i], Results.chord[t, i], Results.beta[t, i]] = [mu,chord, beta]


        a, ap = (0.2,0.02)
        delta = 1e-6
        N_iter_max = 1000
        residual = 1e8
        UR_factor = 0.25
        for ite in range(N_iter_max):
            u_tan, u_nor, u_rel, phi = RelativeVelocities(Rotor, a, ap, r)
            alpha = phi - (beta + Rotor.theta)*np.pi/180
            cl, cd = AirfoilCoefficients(Rotor, alpha)
            lift, drag, f_tan, f_nor = Forces(Rotor, chord, phi, u_rel, cl, cd)

            CT = f_nor*Rotor.n_blades/(0.5*Rotor.rho*Rotor.wind_speed**2*2*m.pi*r)

            if t==0:
                a_new = NewInductionFactor(CT)
            else:
                a_new = Pitt_Peters(CT, Results.a[t-1, i]*Results.f[t-1, i], r, dt, Rotor)

            f,f_tip,f_root = PrandtlTipCorrection(Rotor, mu, a_new)

            a_new = a_new/f
            a = (1-UR_factor)*a + UR_factor*a_new

            if a>0.95 or m.isnan(a_new):
                a=0.95

            ap_new = f_tan*Rotor.n_blades/(2*Rotor.rho*2*m.pi*r*Rotor.wind_speed**2*(1-a)*Rotor.TSR*mu*f)
            ap = (1-UR_factor)*ap + UR_factor*ap_new
            if np.abs(a_new-a) < delta and np.abs(ap_new-ap) < delta:
                break


        #Calculate circulation
        Results.circulation[t, i] = lift/(Rotor.rho*u_rel)


        #Calculate local torque coefficient
        Results.local_CQ[t,i] = f_tan*mu*Rotor.radius*Rotor.n_blades/(0.5*Rotor.rho*Rotor.wind_speed**2*Rotor.radius**2*mu)

        #Store all the results
        [Results.a[t,i], Results.ap[t,i], Results.phi[t,i], Results.alpha[t,i], Results.cl[t,i], Results.cd[t,i],
         Results.f_nor[t,i], Results.f_tan[t,i], Results.f[t,i], Results.f_tip[t,i], Results.f_root[t,i], Results.ite[t,i],
         Results.local_CT[t,i]] = \
        [a, ap, phi*180/np.pi, alpha*180/np.pi, cl, cd, f_nor, f_tan, f, f_tip, f_root, ite, CT]


# #Integrate forces to get total CP, CT, and CQ
Results.Integrate(Rotor_lst, len(time))

# # #Calculate the global axial induction factor
# Results.a_global = NewInductionFactor(self.Results.CT, self.Rotor.yaw)

x = 6  # Want figures to be A6
plt.rc('figure', figsize=[46.82 * .5**(.5 * x), 33.11 * .5**(.5 * x)]   )
#plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plot_TSR(Results,Rotor_lst, [0,1])

plt.figure()
plt.plot(time, Results.CT)
plt.grid()
plt.ylabel('CT [-]')
plt.xlabel('Time [s]')



#%%



