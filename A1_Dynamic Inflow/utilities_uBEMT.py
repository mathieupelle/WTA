# -*- coding: utf-8 -*-
"""
Created on Tue Jun 15 18:54:13 2021

@author: Dell
"""

import numpy as np
#import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import least_squares,fsolve
import math as m
from time import sleep
import warnings


class Rotor:
  def __init__(self, N_radial_sections = 30, Spacing_method = 'lin' ,Optimized_geometry = None):
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
    self.theta = -2 #Pitch angle [deg]
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

    #If there is Optimized geometry, take the chord, pitch, and twist from there
    if Optimized_geometry:
        self.theta = Optimized_geometry.theta*180/np.pi
        self.chord = Optimized_geometry.chord
        self.beta = Optimized_geometry.beta*180/np.pi

    self.SetOperationalData(wind_speed=10, TSR=8, yaw=0) #Assign default values to operational conditions

  def SetOperationalData(self,wind_speed,TSR,yaw,rho=1.225):
    """
      Operational data associated to the rotor

      Parameters
      ----------
      wind_speed : Float [m/s]
      TSR : Float [-]
      yaw : Float [deg]
      rho : Float [kg/m3], optional
           The default is 1.225.
    """
    self.wind_speed = wind_speed
    self.TSR = TSR
    self.yaw = yaw*np.pi/180 #Input yaw should be in degrees!
    self.omega = wind_speed*TSR/self.radius
    self.rho = rho


class Results: #Create the variables to store the results from BEMT
    def __init__(self,N_radial,N_azimuth,N_time):
        self.N_time = N_time
        self.a,self.ap,self.phi,self.alpha,self.cl,self.cd,self.f_nor,self.f_tan,self.f,self.f_tip,self.f_root,self.ite,self.chord,self.beta,self.mu,self.circulation,self.enthalpy_3,self.local_CT,self.local_CQ =  np.zeros((19,N_radial-1,N_azimuth-1,N_time))
        self.azimuth = np.zeros(N_azimuth-1)
        self.CT, self.CP, self.CQ = np.zeros((3, self.N_time))

    def Integrate(self,Rotor,t):
        """


        Parameters
        ----------
        Rotor : Rotor class coresponding to the current time step
        t : Index of the current time step

        Returns
        -------
        None.

        """
        #Calculate global CT
        #Generate differential radius array (to take into account cosinusoidal spacing)
        d_r = np.zeros(len(Rotor.mu)-1)
        for i in range(len(Rotor.mu)-1):
            d_r[i] = (Rotor.mu[i+1]-Rotor.mu[i])*Rotor.radius

        d_azi = 2*np.pi/np.size(self.a,1)

        f_nor = self.f_nor[:,:,t]
        f_tan = self.f_tan[:,:,t]

        self.CT[t] = np.sum(f_nor.transpose()*Rotor.n_blades*d_azi/2/np.pi*d_r)/(0.5*Rotor.rho*Rotor.wind_speed**2*np.pi*Rotor.radius**2)

        #Global CP
        dTorque = f_tan[:,0]*d_r*self.mu[:,0,0]*Rotor.radius*d_azi/2/np.pi
        self.CP[t] = np.sum(dTorque*Rotor.n_blades*Rotor.omega)/(0.5*Rotor.rho*Rotor.wind_speed**3*np.pi*Rotor.radius**2)

        #Global CQ
        self.CQ[t] = np.sum(dTorque*Rotor.n_blades)/(0.5*Rotor.rho*Rotor.wind_speed**2*np.pi*Rotor.radius**3)


class BEMT:
    def __init__(self,Rotor):
        self.Rotor = Rotor

    def RelativeVelocities(self,a,ap,mu,azimuth=0):
        u_tan = self.Rotor.omega*self.Rotor.radius*mu*(1+ap) + self.Rotor.wind_speed*np.sin(self.Rotor.yaw)*np.sin(azimuth)

        #Glauert correction for the normal velocity
        psi = (0.6*a+1)*self.Rotor.yaw
        K = 2*np.tan(psi/2)
        u_nor = self.Rotor.wind_speed*(np.cos(self.Rotor.yaw)-a*(1+K*mu*np.sin(self.Rotor.yaw)))

        #Total relative velocity and flow angle
        u_rel = np.sqrt(u_tan**2 + u_nor**2)
        phi = np.arctan2(u_nor,u_tan)

        return u_tan,u_nor,u_rel,phi


    def AirfoilCoefficients(self,alpha):
        cl = np.interp(alpha*180/np.pi,np.array(self.Rotor.polars['alpha']),np.array(self.Rotor.polars['Cl']))
        cd = np.interp(alpha*180/np.pi,np.array(self.Rotor.polars['alpha']),np.array(self.Rotor.polars['Cd']))

        return cl,cd

    def Forces(self,chord,phi,u_rel,cl,cd):
        lift = 0.5*self.Rotor.rho*u_rel**2*chord*cl
        drag = 0.5*self.Rotor.rho*u_rel**2*chord*cd

        f_tan = lift*np.sin(phi) - drag*np.cos(phi)
        f_nor = lift*np.cos(phi) + drag*np.sin(phi)

        return lift,drag,f_tan,f_nor

    def NewInductionFactor(self,CT,a=None):
        """
        Calculate the new axial induction factor for a given CT and yaw angle.

        Parameters
        ----------
        CT : Value of the thrust coefficient we want to evaluate.
        a : Previous axial induction factor. Used to solve the quadratic equation for the yawed case.
            Default is None.
            If no input is provided, the a_new returned will be the smallest one from the solution of the quadratic equation.
            This is done at this way because we are not using a correction for heavily loaded rotor for the yawed case.

        Returns
        -------
        a_new : Value of the axial induction factor calculated

        """

        if self.Rotor.yaw == 0:
            CT_1 = 1.816 #Constants for the empirical calculation
            CT_2 = 2*np.sqrt(CT_1) - CT_1
            if CT < CT_2:
                a_new = 0.5 - np.sqrt(1-CT)/2
            else:
                a_new = 1 + (CT-CT_1)/(4*np.sqrt(CT_1)-4)

        #Yawed case
        else:
            if a:
                a_new = CT/(4*np.sqrt(1-a*(2*np.cos(self.Rotor.yaw)-a))) #If we have the value from previous iteration use it
            else: #Otherwise, solve numerically
                func = lambda a : 4*a*np.sqrt(1-a*(2*np.cos(self.Rotor.yaw)-a)) - CT
                a_guess = 1/3
                a_new = fsolve(func,a_guess)

        return a_new


    def PrandtlTipCorrection(self,mu,a_new):
        mu_root = self.Rotor.mu[0]
        #Tip correction
        exp = np.exp(-self.Rotor.n_blades/2 * ((1-mu)/mu) * np.sqrt(1+self.Rotor.TSR**2*mu**2/(1-a_new)**2))
        f_tip = 2/np.pi * np.arccos(exp)
        #Root correction
        exp = np.exp(-self.Rotor.n_blades/2 * ((mu-mu_root)/mu) * np.sqrt(1+self.Rotor.TSR**2*mu**2/(1-a_new)**2))
        f_root = 2/np.pi * np.arccos(exp)
        #Combined correction
        f = f_tip*f_root

        if f < 1e-4 or m.isnan(f):
            f = 1e-4

        return f,f_tip,f_root


    def Solver(self, time = [0], conditions = [], DI_Model = "Steady" ,Prandtl_correction = True,N_iter_max = 1000,delta=1e-6,):
        warnings.simplefilter('ignore') #Ignore error messages (division by 0 at the innermost sections with very high nº of points)


        if self.Rotor.yaw == 0:
            N_azimuth = 2
        else:
            N_azimuth = self.Rotor.N_azimuth

        self.Results = Results(self.Rotor.N_radial,N_azimuth,len(time))


        for t_idx in range(len(time)):

            # Set the operational conditions of this time step
            # If steady calculation, takes the values contained in "Rotor"
            if len(time) > 1:
                self.Rotor.SetOperationalData(wind_speed = conditions['wind_speed'][t_idx], TSR = self.Rotor.TSR, yaw = conditions['yaw_angle'][t_idx])
                self.Rotor.theta = conditions['pitch_angle'][t_idx]

                dt = time[1] - time[0] #Get dt

            for i in range(self.Rotor.N_radial-1): #Loop of each blade section
                #Calculate blade parameters in this section
                mu = (self.Rotor.mu[i]+self.Rotor.mu[i+1])/2
                chord = np.interp(mu,self.Rotor.mu,self.Rotor.chord)
                beta = np.interp(mu,self.Rotor.mu,self.Rotor.beta)

                #Store them in the results class as well
                [self.Results.mu[i],self.Results.chord[i],self.Results.beta[i]] = [mu,chord,beta]


                for j in range(N_azimuth-1):
                        azimuth = (self.Rotor.azimuth[j]+self.Rotor.azimuth[j+1])/2
                        self.Results.azimuth[j] = azimuth

                        a,ap = (0.2,0.02) #Initialize induction factors
                        for ite in range(N_iter_max):
                            #Ve m,.-locities and angles
                            [u_tan,u_nor,u_rel,phi] = self.RelativeVelocities(a,ap,mu,azimuth)
                            alpha = phi - (beta + self.Rotor.theta)*np.pi/180

                            #Airfoil forces
                            [cl,cd] = self.AirfoilCoefficients(alpha)
                            [lift,drag,f_tan,f_nor] = self.Forces(chord,phi,u_rel,cl,cd)

                            #Thrust coefficient
                            CT = f_nor*self.Rotor.n_blades/(0.5*self.Rotor.rho*self.Rotor.wind_speed**2*2*np.pi*mu*self.Rotor.radius)

                            #Get new value of axial induction factor
                            if t_idx == 0 or DI_Model == "Steady":
                                 a_new = self.NewInductionFactor(CT,a)
                            else:
                                 if DI_Model == "PP":
                                     a_new = self.Pitt_Peters(CT, self.Results.a[i,j,t_idx-1]*self.Results.f[i,j,t_idx-1], mu*self.Rotor.radius, dt, self.Rotor.wind_speed)
                                 elif DI_Model == "LM":
                                     CT_qs = self.getCT_fromPitchAngle(self.Rotor.theta, self.Rotor.TSR) #Calculate quasi-steady CT based on CT_pitch contours
                                     a_new = self.Larsen_Madsen(self.Results.a[i,j,t_idx-1]*self.Results.f[i,j,t_idx-1], self.Results.local_CT[i,j,t_idx-1], mu*self.Rotor.radius, dt, self.Rotor.wind_speed)
                                 elif DI_Model == "O":
                                     if t_idx==1:
                                         v_int = -self.NewInductionFactor(self.Results.local_CT[i,j,t_idx-1])*self.Rotor.wind_speed

                                     a_new,v_int = self.Oye(a*f, self.Results.local_CT[i,j,t_idx-1], CT, v_int, mu*self.Rotor.radius, dt, self.Rotor.wind_speed)

                                 else:
                                     raise Exception('Its got a C in it. Also you have not selected a model')

                            #Apply the tip and root loss correction factor if wanted
                            if Prandtl_correction:
                                [f,f_tip,f_root] = self.PrandtlTipCorrection(mu,a_new)
                            else:
                                [f,f_tip,f_root] = [1,1,1]
                            a_new = a_new/f

                            #Induction factor for the next iteration
                            a = 0.75*a + 0.25*a_new

                            #Bound a to 0.95
                            if a>0.95 or m.isnan(a_new):
                                a=0.95

                            #Calculate tangential induction
                            ap_new = f_tan*self.Rotor.n_blades/(2*self.Rotor.rho*2*np.pi*mu*self.Rotor.radius*self.Rotor.wind_speed**2*(1-a)*self.Rotor.TSR*mu*f)

                            #Tangential induction for next iteration
                            ap = 0.75*ap + 0.25*ap_new

                            #Check convergency
                            if np.abs(a_new-a) < delta and np.abs(ap_new-ap) < delta:
                                break

                        #Correct values that might be nan due to tip and root effects
                        if m.isnan(f_tan):
                            f_tan = 0
                            f_nor = 0

                        #Calculate circulation
                        self.Results.circulation[i,j,t_idx] = lift/(self.Rotor.rho*u_rel)

                        #Calcualte enthalpy at station 3 of the stream tube
                        self.Results.enthalpy_3[i,j,t_idx] = 1/2*self.Rotor.wind_speed**2*(1-2*a)**2

                        #Calculate local torque coefficient
                        self.Results.local_CQ[i,j,t_idx] = f_tan*mu*self.Rotor.radius*self.Rotor.n_blades/(0.5*self.Rotor.rho*self.Rotor.wind_speed**2*self.Rotor.radius**2*mu)

                        #Store all the results
                        [self.Results.a[i,j,t_idx],self.Results.ap[i,j,t_idx],self.Results.phi[i,j,t_idx],self.Results.alpha[i,j,t_idx],self.Results.cl[i,j,t_idx],
                         self.Results.cd[i,j,t_idx],self.Results.f_nor[i,j,t_idx],self.Results.f_tan[i,j,t_idx],self.Results.f[i,j,t_idx],self.Results.f_tip[i,j,t_idx],self.Results.f_root[i,j,t_idx],self.Results.ite[i,j,t_idx],self.Results.local_CT[i,j,t_idx]] = \
                            [a,ap,phi*180/np.pi,alpha*180/np.pi,cl,cd,f_nor,f_tan,f,f_tip,f_root,ite,CT]

            #Integrate forces to get total CP, CT, and CQ
            self.Results.Integrate(self.Rotor,t_idx)

            #Calculate the global axial induction factor
#            self.Results.a_global = self.NewInductionFactor(self.Results.CT, self.Rotor.yaw)


    def Pitt_Peters(self, CT, a, r, dt, wind_speed):
        """


        Parameters
        ----------
        CT : TYPE
            Current thrust coefficient
        a : TYPE
            Current induction factor
        r : TYPE
            DESCRIPTION.
        dt : TYPE
            DESCRIPTION.
        wind_speed : TYPE
            DESCRIPTION.

        Returns
        -------
        a_new : TYPE
            DESCRIPTION.

        """
        CT_qs = self.NewCT(a)
        dvdt = 3*m.pi*wind_speed**2/16/r*(CT - CT_qs)
        v = -a*wind_speed - dvdt*dt
        a_new = -v/wind_speed
        return a_new

    def Larsen_Madsen(self, a, CT, r, dt, wind_speed):
        vz = -a * wind_speed #calculate induced velocity from a(i)
        tau = 0.5*r/(wind_speed + vz)
        a_qs = self.NewInductionFactor(CT) #calculate a_qs from CT(i+1)
        a_new = a*np.exp(-dt/tau) + a_qs*(1-np.exp(-dt/tau))
        return a_new

    def Oye(self, a1, CT1, CT2, vint1, r, dt, wind_speed):

        # calculate quasi-steady induced velocity
        # vqs1 = a1 * wind_speed
        vqs1 = -self.NewInductionFactor(CT1) * wind_speed

        # calculate current induced velocity
        vz1 = -a1 * wind_speed

        # calculate model time scales
        t1 = (1.1 / (1-1.3*a1)) * (self.Rotor.radius/wind_speed)
        t2 = (0.39-0.26*(r/self.Rotor.radius)**2)*t1 # from Carlos' Jupyter notebook
        # t2 = ((r/self.Rotor.radius)**2)*t1 # from slides

        # next-time-step quasi-steady induction velocity
        vqs2 = -self.NewInductionFactor(CT2) * wind_speed

        # time derivative of intermediate velocity
        dvint_dt = (1/t1) * (vqs1 + 0.6*t1*((vqs2 - vqs1)/dt) - vint1)

        # new intermediate velocity
        vint2 = vint1 + dvint_dt*dt

        # time derivative of induced velocity
        dvz_dt = (1/t2) * (((vint1 + vint2)/2) - vz1)

        # calculate new induced velocity
        vz2 = vz1 + dvz_dt*dt

        # calculate new induction factor
        a2 = -vz2 / wind_speed
        #print(vint1,vz2,CT1,CT2,t1,t2)

        return a2, vint2


    def NewCT(self,a, glauert=True):
        CT = 4*a*(1-a)
        if glauert:
            CT_1 = 1.816
            a_1 = 1-np.sqrt(CT_1)/2
            if a>a_1:
                CT = CT_1-4*(np.sqrt(CT_1)-1)*(1-a)
        return CT


    def CpLambda(self,TSR_list,theta_list):
        """
        Generate Cp-Theta-Lambda contours

        Parameters
        ----------
        TSR_list : Array of TSR values to be analysed.
        theta_list : Array of pitch values to be analysed

        Returns
        -------
        Cp_lambda : Dictionary containing:
            -TSR_list
            -theta_list
            -CP values
            -CT values
            -CQ values
        """

        #Prepare variables to store results
        [CP,CT,CQ] = np.zeros((3,len(TSR_list),len(theta_list)))

        #Store the original pitch angle
        theta_org = self.Rotor.theta

        i=1

        for TSR in TSR_list:
            for theta in theta_list:

                #Assign TSR and theta
                self.Rotor.SetOperationalData(10,TSR,yaw=0)
                self.Rotor.theta = theta

                #Solve BEMT
                self.Solver()

                #Output a status message
                if np.remainder(i,10) == 0:
                    print('Cp-TSR-Theta contours: Point',i,'out of', CP.size,'calculated')
                i = i+1

                #Store the results
                CT[TSR_list.index(TSR),theta_list.index(theta)] = self.Results.CT
                CP[TSR_list.index(TSR),theta_list.index(theta)] = self.Results.CP
                CQ[TSR_list.index(TSR),theta_list.index(theta)] = self.Results.CQ



        #Store all the results in a dictionary for output
        Cp_lambda = {'TSR': TSR_list,
                     'theta': theta_list,
                     'CP':CP,
                     'CT':CT,
                     'CQ':CQ}

        #Also store them in the object for the pitch determination from CT
        self.Cp_lambda = Cp_lambda

        #Assign back the original pitch angle
        self.Rotor.theta = theta_org

        return Cp_lambda

    def getCT_fromPitchAngle(self,theta,TSR):
        #Firstly, check if the Cp-pitch-lambda contours exist
        if hasattr(self, 'Cp_lambda'):
            pass
        else:
            raise Exception("No Cp_lambda variable found. Please run BEMT.CpLambda(TSR_list,theta_list) to generate it so I can interpolate. Thanks")

        #Unpack the necessary variables for readibility
        CT_mat = self.Cp_lambda['CT']
        theta_lst = self.Cp_lambda['theta']
        TSR_lst = self.Cp_lambda['TSR']

        #Interpolate firstly across TSR
        CT_lst = [] #Initialise array of CT vs pitch
        for i in range(len(theta_lst)):
            CT_lst.append(np.interp(TSR,TSR_lst,CT_mat[:,i]))

        #Interpolate the pitch value fiven the desired CT and return it
        return np.interp(theta,theta_lst,CT_lst)



    def getPitchAngle_fromCT(self,CT,TSR):
        """
        Interpolate the pitch angle based on a desired thrust coefficient. It uses the CP-Theta-TSR contours generated before.

        Parameters
        ----------
        CT : Float
            Desired CT which pitch we want to know.
        TSR : Float
            Operational tip speed ratio of the evaluation point.

        Raises
        ------
        Exception
            It needs to have the atribute Cp_lambda generated beforehand.

        Returns
        -------
        Theta: Float [DEGREES]
            Interpolated pitch angle corresponding to the required CT.

        """
        #Firstly, check if the Cp-pitch-lambda contours exist
        if hasattr(self, 'Cp_lambda'):
            pass
        else:
            raise Exception("No Cp_lambda variable found. Please run BEMT.CpLambda(TSR_list,theta_list) to generate it so I can interpolate. Thanks")

        #Unpack the necessary variables for readibility
        CT_mat = self.Cp_lambda['CT']
        theta_lst = self.Cp_lambda['theta']
        TSR_lst = self.Cp_lambda['TSR']

        #Interpolate firstly across TSR
        CT_lst = [] #Initialise array of CT vs pitch
        for i in range(len(theta_lst)):
            CT_lst.append(np.interp(TSR,TSR_lst,CT_mat[:,i]))

        #Sort CT and theta arrays for the second interpolation
        CT_arr = np.array(CT_lst)
        theta_arr = np.array(theta_lst)
        idxs = CT_arr.argsort() #Get indices of sorted array
        CT_arr = CT_arr[idxs]
        theta_arr = theta_arr[idxs]

        #Interpolate the pitch value fiven the desired CT and return it
        return np.interp(CT,CT_arr,theta_arr)


class Optimizer:
    def __init__(self, Rotor_original, a, TSR):
        self.a = a
        self.R = Rotor_original.radius
        self.TSR = TSR
        self.B = Rotor_original.n_blades
        self.mu = Rotor_original.mu

        #Calculate optimal Cl and E
        Alpha = Rotor_original.polars['alpha']
        Cl = Rotor_original.polars['Cl']
        Cd = Rotor_original.polars['Cd']

        #Select the point with the maximum efficiency
        self.E = max(Cl/Cd)
        self.cl = Cl[np.argmax(Cl/Cd)]
        self.cd = Cd[np.argmax(Cl/Cd)]
        self.aoa = Alpha[np.argmax(Cl/Cd)]

        #Execute the optimization for chord and twist
        self.ChordOpt()
        self.TwistOpt()


    def residuals(self,x):

        c,ap = x #Unpack the input

        #Flow angle
        phi = m.atan((1-self.a)*self.R/((1+ap)*self.r*self.TSR))

        #Tip loss
        exp = np.exp(-self.B/2 * ((1-self.r/self.R)/(self.r/self.R)) * np.sqrt(1+self.TSR**2*(self.r/self.R)**2/(1-self.a)**2))
        f_tip = 2/np.pi * np.arccos(exp)
        #Root correction
        exp = np.exp(-self.B/2 * ((self.r/self.R-0.2)/(self.r/self.R)) * np.sqrt(1+self.TSR**2*(self.r/self.R)**2/(1-self.a)**2))
        f_root = 2/np.pi * np.arccos(exp)
        ##Combined correction
        F = f_tip*f_root

        #Force coefficients
        Cy = self.cl * np.cos(phi) + self.cd*np.sin(phi)
        Cx = self.cl * np.sin(phi) - self.cd*np.cos(phi)

        #Solidity
        sigma = c*self.B/(2*np.pi*self.r)

        #Get residual c and ap
        res_c = 4*np.pi*self.r*m.sin(phi)**2*F*2*self.a/(Cy*self.B*(1-self.a)) - c
        res_ap = 1/(4*F*np.sin(phi)*np.cos(phi)/(sigma*Cx)-1) - ap

        return res_c,res_ap


    def ChordOpt(self):
        [self.chord,self.ap] = np.zeros((2,len(self.mu)))
        for i in range(len(self.mu)):
            self.r = self.mu[i]*self.R #Radial position
            x0 = [3,0.001] #Initial guess
            #Lower and upper bounds
            if self.mu[i]<0.5:
                bounds = ((1.5,0),(7,1)) #Minimum 1.5m at the root for structural reasons
            else:
                bounds = ((0.3,0),(7,1)) #Minimum 0.3 at the tip to avoid c=0
            results = least_squares(self.residuals,x0,bounds=bounds,verbose=0) #Calculate with the least-squares method the chord and a'
            self.chord[i],self.ap[i] = results.x

    def TwistOpt(self):
        self.beta = np.zeros((len(self.mu)))
        for i in range(len(self.mu)):
            r = self.mu[i]*self.R

            #Calculate flow angle at each section
            phi = m.atan((1-self.a)*self.R/((1+self.ap[i])*r*self.TSR))

            #Set the twist such that the optimal angle of attack is seen at each section
            self.beta[i] = phi - self.aoa*np.pi/180

        #When the twist distribution has been calcualted, set the pitch such that the twist at the last section is zero
        self.theta = self.beta[-1]
        self.beta = self.beta - self.theta


def MeshSensitivity(N_array,Spacing):
    [CT,execution_time] = np.zeros((2,len(N_array))) #Initialize an array to store thurst and execution time

    #Loop to each number of points
    for idx in range(len(N_array)):
        print('Solving rotor with',N_array[idx],'radial points. Case',idx+1,'out of',len(N_array)) #Status message
        #Create the rotor with the corresponding nº of elements and spacing method
        Rotor_spacing = Rotor(N_radial_sections = N_array[idx], Spacing_method = Spacing)
        #Start the BEMT class with them
        BEMT_spacing = BEMT(Rotor_spacing)
        #Solve the linear-spaced
        tic = time.perf_counter() #Initialize a timer for the execution time analysis
        BEMT_spacing.Solver() #Solve the rotor
        execution_time[idx] = time.perf_counter() - tic #Stop the count!
        #Solve the cosinusioidal-spaced
        CT[idx] = BEMT_spacing.Results.CT #Save thrust coefficient

    #Calculate relative error for each position
    err = abs(CT-CT[-1])/CT[-1]

    #Get the number of mesh points that gives us 0.001 relative error (0.1%)
    N_chosen = np.interp(x = 1e-3,xp = np.flip(err),fp = np.flip(N_array))

    return CT,err,N_chosen,execution_time


save = False #Flag for saving plots

def plot_optimized_geometry(Rotor_org,Rotor_opt,Results_org,Results_opt,Cp_lambda_org,Cp_lambda_opt):

    #Set default stuff
    x = 6  # Want figures to be A6
    plt.rc('figure', figsize=[46.82 * .5**(.5 * x), 33.11 * .5**(.5 * x)]   )
    #plt.rc('text', usetex=True)
    plt.rc('font', family='serif')

    #Compare blade geometries
    fig = plt.figure('Chord distribution')
    plt.plot(Rotor_org.mu,Rotor_org.chord)
    plt.plot(Rotor_opt.mu,Rotor_opt.chord)
    plt.legend(['Original','Optimized'])
    plt.xlabel('Radius r/R [-]')
    plt.ylabel('Chord [m]')
    plt.grid()

    if save==True:
        plt.savefig('figures/Optimization/Chord.pdf')

    fig = plt.figure('Twist distribution')
    plt.plot(Rotor_org.mu,Rotor_org.beta)
    plt.plot(Rotor_opt.mu,Rotor_opt.beta)
    plt.legend(['Original','Optimized'])
    plt.xlabel('Radius r/R [-]')
    plt.ylabel('Twist [deg]')
    plt.grid()

    if save==True:
        plt.savefig('figures/Optimization/Twist.pdf')

    #Plot CP-lambda-theta contours
    fig = plt.figure('CP-lambda-theta')
    CS = plt.contour(Cp_lambda_org['TSR'],Cp_lambda_org['theta'],Cp_lambda_org['CP'].transpose(),cmap='jet',levels=10)
    plt.clabel(CS, inline=1, fontsize=10)
    plt.xlabel('Tip speed ratio [-]')
    plt.ylabel('Pitch angle [deg]')
    plt.title('Power coefficient CP [-] (Original design)')
    if save==True:
        plt.savefig('figures/Optimization/Cp_Lambda_Theta_org.pdf')

    fig = plt.figure('CP-lambda-theta optimized')
    CS = plt.contour(Cp_lambda_opt['TSR'],Cp_lambda_opt['theta'],Cp_lambda_opt['CP'].transpose(),cmap='jet',levels=10)
    plt.clabel(CS, inline=1, fontsize=10)
    plt.xlabel('Tip speed ratio [-]')
    plt.ylabel('Pitch angle [deg]')
    plt.title('Power coefficient CP [-] (Optimized design)')
    plt.plot(8,Rotor_opt.theta,'x',color='black')
    if save==True:
        plt.savefig('figures/Optimization/Cp_Lambda_Theta_opt.pdf')

    #Plot CT-lambda-theta contours
    fig = plt.figure('CT-lambda-theta')
    CS = plt.contour(Cp_lambda_org['TSR'],Cp_lambda_org['theta'],Cp_lambda_org['CT'].transpose(),cmap='jet',levels=10)
    plt.clabel(CS, inline=1, fontsize=10)
    plt.xlabel('Tip speed ratio [-]')
    plt.ylabel('Pitch angle [deg]')
    plt.title('Thrust coefficient CT [-] (Original design)')
    if save==True:
        plt.savefig('figures/Optimization/Ct_Lambda_Theta_org.pdf')

    fig = plt.figure('CT-lambda-theta optimized')
    CS = plt.contour(Cp_lambda_opt['TSR'],Cp_lambda_opt['theta'],Cp_lambda_opt['CT'].transpose(),cmap='jet',levels=10)
    plt.clabel(CS, inline=1, fontsize=10)
    plt.xlabel('Tip speed ratio [-]')
    plt.ylabel('Pitch angle [deg]')
    plt.title('Thrust coefficient CT [-] (Optimized design)')
    plt.plot(8,Rotor_opt.theta,'x',color='black')
    if save==True:
        plt.savefig('figures/Optimization/Ct_Lambda_Theta_opt.pdf')

    #CQ-lambda-theta
    fig= plt.figure('CQ-lambda-theta')
    CS = plt.contour(Cp_lambda_org['TSR'],Cp_lambda_org['theta'],Cp_lambda_org['CQ'].transpose(),cmap='jet',levels=10)
    plt.clabel(CS, inline=1, fontsize=10)
    plt.xlabel('Tip speed ratio [-]')
    plt.ylabel('Pitch angle [deg]')
    plt.title('Torque coefficient CQ [-] (Original design)')
    if save==True:
        plt.savefig('figures/Optimization/Cq_Lambda_Theta_org.pdf')

    fig = plt.figure('CQ-lambda-theta optimized')
    CS = plt.contour(Cp_lambda_opt['TSR'],Cp_lambda_opt['theta'],Cp_lambda_opt['CQ'].transpose(),cmap='jet',levels=10)
    plt.clabel(CS, inline=1, fontsize=10)
    plt.xlabel('Tip speed ratio [-]')
    plt.ylabel('Pitch angle [deg]')
    plt.title('Torque coefficient CQ [-] (Optimized design)')
    plt.plot(8,Rotor_opt.theta,'x',color='black')
    if save==True:
        plt.savefig('figures/Optimization/Cq_Lambda_Theta_opt.pdf')

    #Compare original vs optimal design
    rho = Rotor_org.rho
    wind_speed = Rotor_org.wind_speed
    radius = Rotor_org.radius
    n_blades = Rotor_org.n_blades
    omega = Rotor_org.omega

    var=['alpha','phi','a','ap','f_tan','f_nor','circulation','local_CQ','f','cl','cd']
    labels=[r'$\alpha$ [deg]','$\phi$ [deg]', 'a [-]','$a^,[-]$', '$C_t$ [-]', '$C_n$ [-]','$\Gamma$ [-]','$C_q [-]$','Prantl\'s tip loss factor','$C_l$','$C_d$']
    for i in range(len(var)):
        plt.figure()
        plt.grid()
        plt.xlabel(r'Radius $\frac{r}{R}$ [-]')
        plt.ylabel(labels[i])

        dic=Results_org['TSR8_yaw0']
        dic_opt = Results_opt['TSR8_yaw0']
        if var[i]=='f_tan' or var[i]=='f_nor':
            Z=getattr(dic, str(var[i]))/(0.5*rho*wind_speed**2*radius)
            Z_opt=getattr(dic_opt, str(var[i]))/(0.5*rho*wind_speed**2*radius)
        elif var[i]=='circulation':
            Z=getattr(dic, str(var[i]))/((np.pi*wind_speed**2/(n_blades*omega)))
            Z_opt=getattr(dic_opt, str(var[i]))/((np.pi*wind_speed**2/(n_blades*omega)))
        else:
            Z=getattr(dic, str(var[i]))
            Z_opt=getattr(dic_opt, str(var[i]))

        plt.plot(dic.mu,Z,label='Original design')
        plt.plot(dic_opt.mu,Z_opt,label='Optimized design')
        plt.legend()
        if save:
            plt.savefig('figures/Optimization/'+str(var[i])+'.pdf')

    #Plot torque differential
    plt.figure()
    plt.grid()
    plt.xlabel(r'Radius $\frac{r}{R}$ [-]')
    plt.ylabel('$\delta Q / (0.5 \\rho U^2 R^2)$')
    Z=getattr(dic, 'f_tan')*dic.mu*radius*n_blades/((0.5*rho*wind_speed**2*radius**2))
    Z_opt=getattr(dic_opt, 'f_tan')*dic.mu*radius*n_blades/((0.5*rho*wind_speed**2*radius**2))

    plt.plot(dic.mu,Z,label='Original design')
    plt.plot(dic_opt.mu,Z_opt,label='Optimized design')
    plt.legend()
    if save:
        plt.savefig('figures/Optimization/Torque_diff.pdf')



def plot_mesh_sensitivity(N_array,CT_lin,CT_cos,N_chosen_lin,N_chosen_cos,execution_time_lin,execution_time_cos,err_lin,err_cos):

    fig,ax = plt.subplots(figsize=(10,5))
    ax.plot(N_array,CT_lin,label='CT (lin)')
    ax.plot(N_array,CT_cos,'--',color='tab:blue',label='CT (cos)')
    ax.set_ylabel('Thrust coefficient $CT$ [-]')
    ax.grid()
    ax2=ax.twinx()
    ax2.plot(N_array,execution_time_lin,color='tab:orange',label='time (lin)')
    ax2.plot(N_array,execution_time_cos,'--',color='tab:orange',label='time (cos)')
    ax2.plot([N_chosen_lin,N_chosen_lin],[execution_time_lin.min(),execution_time_cos.max()],'-',color='black',label='N (lin)')
    ax2.plot([N_chosen_cos,N_chosen_cos],[execution_time_lin.min(),execution_time_cos.max()],'--',color='black',label='N (cos)')
    ax2.set_ylabel('BEMT execution time [s]')
    ax.set_xlabel('Number of radial elements $N$')
    ax.legend(bbox_to_anchor=(1.07,0.7), loc="upper left")
    ax2.legend(bbox_to_anchor=(1.07,0.5), loc="upper left")
    fig.tight_layout()
    fig.subplots_adjust(right=0.8)
    if save:
        plt.savefig('figures/MeshSensitivity/CT_time.pdf')


    fig = plt.figure()
    plt.loglog(N_array[:-1],err_lin[:-1],'-',label='linear spacing')
    plt.loglog(N_array[:-1],err_cos[:-1],'-',label='cosinusoidal spacing')
    plt.grid(which='minor')
    plt.grid()
    plt.plot([N_array.min(),N_array.max()],[1e-3,1e-3],'--',color='black')
    plt.xlabel('Number of radial elements $N$')
    plt.ylabel('Relative error')
    plt.legend()
    if save:
        plt.savefig('figures/MeshSensitivity/Convergence.pdf')

