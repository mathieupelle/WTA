# -*- coding: utf-8 -*-
"""
Created on Tue Jun 15 19:28:40 2021

@author: Dell
"""
import numpy as np
import matplotlib.pyplot as plt
from utilities_DI import Rotor,BEMT,Optimizer,MeshSensitivity,plot_optimized_geometry,plot_mesh_sensitivity
    

#TESTING
Test = Rotor(N_radial_sections = 40)
Test.SetOperationalData(10, 10, 0)

cond = {'wind_speed': [10, 10],
        'pitch_angle': [-2,0],
        'yaw_angle': [0,0]}

time_array = [0,.1]
Sol = BEMT(Test)
Sol.Solver(time=time_array,conditions=cond,DI_Model="PP")

#Testing the CP-LAMBDA contours
Cont = Sol.CpLambda([8,9,10,11], [-2,-1,0,1,2])

pitch = Sol.getPitchAngle_fromCT(0.7, 8)
    

    
plt.figure()
plt.plot(Sol.Results.a[:,0,0])
plt.plot(Sol.Results.a[:,0,1])
