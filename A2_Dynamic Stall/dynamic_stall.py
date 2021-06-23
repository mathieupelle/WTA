# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

#   Import python packages to be used
import matplotlib.pyplot as plt
# import matpLotLib.coLors as coLors
#_from mpL_tooLkits.mpLot3d import Axe53D
import numpy as np
#_fr0m matpLotLib.transfbrms import Bbox
# import matpLotLib

#------------------------

#   Define conversion between reduced time and reduced frequency
def time2semichord(time):
    return 2*Uinf*time/chord

def semichord2time(s):
    return s/2/Uinf*chord

#------------------------

#   Define the case study
# define properties of the system
dt=.1
time =np.arange(0,500,dt)
Uinf=1
# properties Qf the airfoiL
chord=1 # chord of the airfoiL
a_dalpha=2*np.pi # Lift sLope
alpha0=0 # aLpha for which normaL Load is zero in steady_fLow
# pitching motion of the airfoiL
k=.1 # reduced_frequency of the pitching motion
omega=k*2*chord/Uinf #_frequency Qf the piching motion
Amplitude_alpha=10/180*np.pi # ampLitude of the pitching motion
alpha_t0=15/180*np.pi # aLpha at time=9
alpha=Amplitude_alpha*np.sin(omega*time)+alpha_t0 # caLcuLate aLpha
dalpha_dt=np.gradient(alpha,time) # caLcuLate the time derivative of aLpha
# pLunge motion of the airfoiL
k_plg=.0 # reduced frequency of the pLunge motion
omega_plg=k_plg*2*chord/Uinf #_frequency of the pLunge motion
Amplitude_plg=.3 # ampLitude of the pLunge motion
hplg=Amplitude_plg*np.sin(omega_plg*time) #position
dhplg_dt=np.gradient(hplg,time) # pLunge veLocity
# define the array semi-chord time scaLe
sarray = time2semichord(time)

#------------------

#   Unsteady attached flow model
# caLCULate quasi-steady aLpha
alpha0 = 0 # we define the \aLpha_{0}, zero for the case of an uncambered pLate/ainfoiL
alphaqs = alpha + dalpha_dt*(chord/2)/Uinf - dhplg_dt/Uinf
dalphaqs_dt=np.gradient(alphaqs,time) # caLcuLate the time derivative of the quasi-steaddy aLpha
# caLcuLate the coefficient of normaL_force assuming quasi-steady_fLow asuming potentiaL fLow
Cnormal_quasisteady = 2*np.pi*(alphaqs-alpha0)

# we pLot the effective quasi—steady angLe of attack \aLpha_{qs}
# pLot figure
plt.rcParams.update({'font.size': 14}) #, 'figure.dpi':159, 'savefig.dpi':159})
plt.rcParams["font.family"] = "serif" # define font
plt.rcParams["mathtext.fontset"] = "dejavuserif" # define font
cmap = plt.get_cmap('BuGn') # define coLormap
fig,ax = plt.subplots(figsize=[6,6]) # define pointers for the figure and axes
ax.plot(alpha*180/np.pi, alphaqs*180/np.pi,color='black', linewidth=1) # pLot equivaLent quasi-steady angLe of attack
ax.set_xlabel(r'$\alpha (^\circ)$') # set x-LabeL
ax.set_ylabel(r'$\alpha_{qs} (^\circ)$') # set y-LabeL
# add arrows to indicate the direction of the cyCLe
parr1=ax.annotate('', xy=(17.5, 20), xytext=(10,12.5),
                  arrowprops=dict(color='black', shrink=0.05, width=.5, headwidth=3,headlength=4, linewidth=.2))
parr1=ax.annotate('', xy=(10, 7.5), xytext=(17.7,15),
                  arrowprops=dict(color='black', shrink=0.05, width=.5, headwidth=3,headlength=4, linewidth=.2))
plt.grid() # add a grid
ax.set_xlim(0,30) # define Limits of the axis
ax.set_ylim(0,30) # define Limits of the axis
plt.tight_layout() # aLL eLements of_figure inside pLot area
plt.show() # show_figure
filename = 'figures/alpha_quasi_steady' # define name of the figure to be saved
#fig.savefig(filename+'.svg', pad_inches = 0) # save figure
#fig.savefig(filename+'.pdf',pad_inches = 0) # save_figure
#fig.savefig(filename+'.png', pad_inches = 0, dpi=300) # save_figure

# we pLot the quasi—steady normaL—force coefficient
# pLot figure
plt.rcParams.update({'font.size': 14}) #, 'figure.dpi':159, 'savefig.dpi':159})
plt.rcParams["font.family"] = "serif" # define font
plt.rcParams["mathtext.fontset"] = "dejavuserif" # define font
cmap = plt.get_cmap('BuGn') # define coLormap
fig,ax = plt.subplots(figsize=[6,6]) # define pointers for the figure and axes
ax.plot(alpha*180/np.pi, Cnormal_quasisteady,color='black', linewidth=1) # pLot equivaLent quasi-steady angLe of attack
ax.set_xlabel(r'$\alpha (^\circ)$') # set X-LabeL
ax.set_ylabel(r'$Cn_{qs} $') # set y-LabeL
# add arrows to indicate the direction of the cycLe
parr1=ax.annotate('', xy=(17.5, 20/360*4*np.pi**2), xytext=(10,12.5/360*4*np.pi**2),
                  arrowprops=dict(color='black', shrink=0.05, width=.5, headwidth=3,headlength=4, linewidth=.2))
parr1=ax.annotate('', xy=(10, 7.5/360*4*np.pi**2), xytext=(17.7,15/360*4*np.pi**2),
                  arrowprops=dict(color='black', shrink=0.05, width=.5, headwidth=3,headlength=4, linewidth=.2))
plt.grid() # add a grid
ax.set_xlim(0,30) # define Limits of the axis
ax.set_ylim(0,3) # define Limits of the aXis
plt.tight_layout() # aLL eLements of figure inside pLot area
plt.show() # show_figure
filename = 'figures/Cnormal_quasisteady' # define name of the figure to be saved
#fig.savefig(filename+'.svg', pad_inches = 0) # save_figure
#fig.savefig(filename+'.pdf',pad_inches = 0) # save_figure
#fig.savefig(filename+'.png', pad_inches = 0, dpi=300) # save_figure

#-------------------------

#   Calculation of the unsteady normal force

# determining X and Y terms for recursive marching formuLa for approximation of DuhameL's integraL
def duhamel_approx(Xi,Yi,delta_s,delta_alpha,order=2,A1=0.3,A2=0.7,b1=0.14,b2=0.53):
    # A1=9.165,A2=0.335,b1=0.0455,b2=0.3
    # determine the next vaLues ofAX and Y; named Xipl and Yipl
    if order==1:
        Xip1= Xi*np.exp(-b1*delta_s)+A1*delta_alpha
        Yip1= Yi*np.exp(-b2*delta_s)+A2*delta_alpha
    elif order==2:
        Xip1= Xi*np.exp(-b1*delta_s)+A1*delta_alpha*np.exp(-b1*delta_s/2)
        Yip1= Yi*np.exp(-b2*delta_s)+A2*delta_alpha*np.exp(-b2*delta_s/2)
    else:
        Xip1= Xi*np.exp(-b1*delta_s)+A1*delta_alpha*((1+4*np.exp(-b1*delta_s/2)+np.exp(-b1*delta_s))/6)
        Yip1= Yi*np.exp(-b2*delta_s)+A2*delta_alpha*((1+4*np.exp(-b2*delta_s/2)+np.exp(-b2*delta_s))/6)
    return Xip1,Yip1
# define_function for circuLatory force, potentiaL_fLow
def circulatory_normal_force(a_dalpha,alpha_equivalent,alpha0):
    return a_dalpha*(alpha_equivalent-alpha0)

# define arrays for X;Y and aLpha_equivaLent
Xarray=np.zeros(np.shape(time))
Yarray=np.zeros(np.shape(time))
# define the array of aLpha_eauivaLent
alpha_equivalent=np.zeros(np.shape(time))
alpha_equivalent[0]=alphaqs[0]
# march soLution in time for aLpha_E
for i,val in enumerate(time[:-1]):
    Xarray[i+1],Yarray[i+1]=duhamel_approx(Xarray[i],Yarray[i],sarray[i+1]-sarray[i],alphaqs[i+1]-alphaqs[i])
alpha_equivalent=alphaqs-Xarray-Yarray
#print("alpha_equivalent",alpha_equivalent)

# pLot soLutions of test of duhameL_appr0x
# pLot_figure
plt.rcParams.update({'font.size': 14}) #, 'figure.dpi':159, 'savefig.dpi':159})
plt.rcParams["font.family"] = "serif" # define fbnt
plt.rcParams["mathtext.fontset"] = "dejavuserif" # define fbnt
cmap = plt.get_cmap('BuGn') # define coLormap
fig,ax = plt.subplots(figsize=[6,6]) # define pointers fbr figure and axes
#we wiLL onLy pLot the Last cycLe
Ncycles = np.floor(time[-1]*omega/(2*np.pi)) # determine number of eycLes
n_of_cycle = time*omega/(2*np.pi) # caLcuLate the phase Qf the different points Qf the cycLe
i1=np.argmin(np.abs(n_of_cycle-(Ncycles-1))) # index af start of eycLe pLotted
i2=np.argmin(np.abs(n_of_cycle-(Ncycles-.5))) # index of 189 degrees
i3=np.argmin(np.abs(n_of_cycle-(Ncycles))) # index of'369 degrees
# pLot Last cycLe of the simuLation, steady, quasi—steady and unsteady eauivaLent angLe of attack
ax.plot(time2semichord(n_of_cycle[i1:i3]-n_of_cycle[i1]), alpha[i1:i3]*180/np.pi,color='blue',linestyle='--', label=r'$\alpha$')
ax.plot(time2semichord(n_of_cycle[i1:i3]-n_of_cycle[i1]), alphaqs[i1:13]*180/np.pi,color='red',linestyle='—.', label=r'$\alpha_{qs}$')
ax.plot(time2semichord(n_of_cycle[i1:i3]-n_of_cycle[i1]), alpha_equivalent[i1:i3]*180/np.pi,color='green',linestyle='-', label=r'$\alpha_{eq}$')
ax.set_xlabel('s semichords') # set x-LabeL
ax.set_ylabel(r'$(6 \circ)$') # set y-LabeL
ax.set_xlim(0,2) # define Limits af the axis
ax.set_ylim(0,30) # define Limits af the axis
ax.grid() # add grid
ax.legend(loc='lower left')
plt.tight_layout() # aLL eLements af_figure inside pLot area
plt.show() # sh0w_figure
filename = 'figures/comparison_alpha_st_qs_circ' # define name of the figure to be saved
#fig.savefig(filename+'.svg', pad_inches = 0) # save figure
#fig.savefig(filename+'.pdf',pad_inches = 0) # save_figure
#fig.savefig(filename+'.png', pad_inches = 0, dpi=360) # save figure



