#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 20 12:11:03 2026

@author: maggie
"""

"""
Transport Model Version 5 --> POS 3 data!!!!

Hypothesis:
Compartment growth in Xenopus cell-free extract alters local cytoplasmic
transport, which changes the phase lag of CDK1 signaling.

Workflow:

Experimental polymerization rate
        
Cluster diameter
        
Cluster volume
        
Cluster growth velocity
        
Local cytoplasmic velocity
        
Transport coefficient
        ↓
Predicte phase lag

Version 5 assumptions:
- Baseline cytoplasmic velocity is normalized to 1.
- Low Reynolds number (Stokes flow).
- No-slip boundary condition at the compartment boundary.
- Pressure gradient represented by a normalized constant.

NOTE: refer to "physical biology of the cell" book when you get stuck with equations
      and print out those plots!!! you'll see where there is a scalling issue faster!!!
"""

#Library Calls ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
import numpy as np
import random
import matplotlib.pyplot as plt
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#Global Variables ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
t = 1200 #mins
dt = 1 #time step, time in between frames
cyto_visc = 1.20e-03 #known visocity of cytoplasm, subject to change, we could get our own data with rheometer
T = 60 #constant period, 60 minutes, typicall of cell-free Xenopus extract 
pressure_grad = 1.0 #normalized for now, will change later if neccesary 
vel_baseline = 1.0 #normalized but will change later
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#Functions ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def volume_calc(radius):
    cluster_volume = (4/3) * np.pi * (radius/2)**3
    return cluster_volume #micron^3

#HAVE TO ADD IN VELCOIDTY CALC BEFORE ULOCAL CALC 

def cluster_vel_calc(d_radius):
    if d_radius == 0:
        cluster_growth_rate = vel_baseline
    else:
        cluster_growth_rate = d_radius / dt
    return cluster_growth_rate

def localvel(velgrowth, dynamic_visc):
    if velgrowth > vel_baseline:
        ulocal = (1/dynamic_visc) * pressure_grad
    else:
        ulocal = vel_baseline

def transport_coe_calc(vclust):
    T_coe = vclust/vel_baseline
    return T_coe

def delta_phase(transcoe):
    d_phase = 1 / transcoe
    return d_phase
    
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#Bins, data collectors ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Min = []
curr_diameter = []
curr_radius = []
change_in_radii = []
curr_volume = []
cluster_vel = []
cyto_vel_at_cluster = []
transport_coeffiecent = []
phase_lag_change = []

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


current_cycle = 0
cycle_time = 0


#Polyrate data, measured in ImageJ ~~~~~~~~~~~~~~~~~~~~~~
#loop updating diameter based on the poly rates ~~~~~~~~~
for i in range(0, t, dt):
    
    current_cycle = i // T
    cycle_time = i % T
    
    
    if current_cycle == 0:
        effective_diameter = 0.008
        effective_radius = effective_diameter / 2
        polyrate = 0
    
    else:
        
        #cycles 1-3   
        if 1 <= current_cycle < 3:
            if cycle_time == 0:
                effective_diameter = 0.008
                effective_radius = effective_diameter/2
                polyrate = 0
            elif 1 <= cycle_time <= 35:
                polyrate = (10.583 * cycle_time - 68.421) * 2.6
                if polyrate > effective_diameter:
                    effective_diameter = polyrate #* random.uniform(.99, 1.01)
                    effective_radius = effective_diameter/2
            elif 35 < cycle_time < 60:
                polyrate = (-5.3734 * cycle_time + 352.23) * 2.6
                decay_effective_diameter = polyrate #* random.uniform(.99, 1.01)
                effective_diameter = decay_effective_diameter
                if effective_diameter < 0.008:
                    effective_diameter = 0.008
                effective_radius = effective_diameter/2
    
        #cycles 3-8
        elif 3 <= current_cycle < 9:
            if cycle_time == 0:
                effective_diameter = 0.008
                effective_radius = effective_diameter/2
                polyrate = 0
            elif 1 <= cycle_time <= 24:
                polyrate = (7.6703 * cycle_time - 22.956) * 2.6
                if polyrate > effective_diameter:
                    effective_diameter = polyrate #* random.uniform(.99, 1.01)
                    effective_radius = effective_diameter/2
            elif 24 < cycle_time < 60:
                polyrate = (0.8247 * cycle_time + 174.46) * 2.6
                if polyrate > effective_diameter:
                    effective_diameter = polyrate #* random.uniform(.99, 1.01)
                    effective_radius = effective_diameter/2
    
        #cyles 9 and past
        else:
            if cycle_time == 0:
                effective_diameter = 0.008
                effective_radius = effective_diameter/2
                polyrate = 0
            elif 1 <= cycle_time <= 38:
                polyrate = ((6.7632 * cycle_time) - 6.5) * 2.6
                if polyrate > effective_diameter:
                    effective_diameter = polyrate #* random.uniform(.99, 1.01)
                    effective_radius = effective_diameter/2
            elif 38 < cycle_time < 60:
                polyrate = (0.0296 * cycle_time + 212.39) * 2.6
                if polyrate > effective_diameter:
                    effective_diameter = polyrate #* random.uniform(.99, 1.01)
                    effective_radius = effective_diameter/2
                
        
# End of Diameter and Polyrate calcs~~~~~~~~~~~~~~~~~~~~~~~
    
    #add data to the bin
    Min.append(i)
    curr_diameter.append(effective_diameter)
    curr_radius.append(effective_radius)
    
    #volume calculation~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    calculated_volume = volume_calc(effective_radius)
    #add data to bin 
    curr_volume.append(calculated_volume)
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    #change in radii~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if cycle_time == 0:
        delta_radii = 0

    else:
        delta_radii = curr_radius[i] - curr_radius[i-1] #something is wrong here:
            #I don't think the change in radii is being found:
                #ex: time 67-69: delta radii = 7..358 --> 13.7579 --> 21.1159
    #add to bin 
    change_in_radii.append(delta_radii)
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    #the velocidty of teh tubulin cluster ~~~~~~~~~~~~~~~~~
    #important for caluclating the change in vel of cyto later
    #refer to the book: "physical biology of the cell" for more info...chpt 12.4
    cluster_growth_velocity = cluster_vel_calc(delta_radii)
    #add data to bin 
    cluster_vel.append(cluster_growth_velocity)
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    #velocity at cluster ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #not sure if this is right....
    cyto_velocity_on_cluster = velofclust(cluster_growth_velocity)
    #add to bin 
    cyto_vel_at_cluster.append(cyto_velocity_on_cluster)
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    #transport coe calc ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    transport_coe = transport_coe_calc(cyto_velocity_on_cluster)
    #add to bin
    transport_coeffiecent.append(transport_coe)
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    #change in phase lag ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    phase_change = delta_phase(transport_coe)
    #add to data bin 
    phase_lag_change.append(phase_change)
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    

    
#END OF CALC BLOCK~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
for j in range(0,t,dt):
    print("MINUTE: ", Min[j], 
          "Radii: ", curr_radius[j],
          "Diameter: ", curr_diameter[j],
          #"Cluster Growth Velocity: ", cluster_vel[j],
          #"Cytoplasm Flow at Cluster: ", cyto_vel_at_cluster[j], 
          "Change in Radii: ", change_in_radii[j]) 
          #"Transport Coeffiecent: ", transport_coeffiecent[j],
          #"Change in Phase Lag/Lead: ", phase_lag_change[j])
    
plt.plot(Min, curr_diameter, label = 'diameter')
plt.legend()  
plt.show()

plt.plot(Min, curr_radius, label = 'radii')
plt.legend()  
plt.show()

plt.plot(Min, change_in_radii, color = 'green', label = 'delta radii') #issue here with first 3 cycles
plt.legend()  
plt.show()

plt.plot(Min, cluster_vel, color = 'm', label = 'cluster growth velocity') #this is the same at above so graph will be same 
plt.legend()  
plt.show()

plt.plot(Min, cyto_vel_at_cluster, color = 'blue', label = 'cytoplasm velocity at cluster') #since our baseline is 1, and delta_radi varies by 100,
                                    #the data isn't changing much, affecting down stream 
plt.legend()                              
plt.show()

plt.plot(Min, transport_coeffiecent, color = 'darkorange', label = 'transport coe' )
plt.legend()  
plt.show()

plt.plot(Min, phase_lag_change, color = 'red', label = 'phase lag')
plt.legend()  
plt.show()