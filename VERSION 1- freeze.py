#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 21 17:09:27 2026

@author: maggie
"""

#Model 1 --> iteration 6 

#Hypothesis: 
    #The compartments formed in cell-free xenopus extract create an enviroment 
    #that favors transportation of cdk1 signals. By treating the extract as a low-reynolds
    #number, the geometry of the compartments directly impacts the relative mitotic phase-lag.
    
#WORK FLOW:
    
    #Diameter --> from polyrate (ImageJ), Experimental (GEOMETRY COMPONENT)
    #Radius --> calculated
    #Growth rate --> calculated (dr/dt), sets boundary condition
    #Occupancy Fraction --> caluclated, from curr_radii / MAX radii (ImageJ)
    #Dynamic Viscoisty --> baseline * e ^-(k*occupation fraction)
    #Effective Diffusion --> baseline_diffusion . visocisity 
    #local flow --> pressure gradient / viscoisty WITH BOUNDARY CONDITION
    #Transport coe --> (ulocal/ubaseline) * (Deff/Do) 
    #Delta Phase Lag --> 1/T
    
#Assumptions 
    #vel_baseline = 1.0
    #pressure_gradient = 1.0
    #Baseline_diffusion = 1.0
    #Assuming that cytoplasm will be treated as a low reynolds number
        #this implies --> no interia 
    #dynamic viscosity calc 
    #occupation fraction calc 
    #k = 1.0, crowding constant 
    #assuming miotic activity is lagging 


#Import Libraries 
import numpy as np
import matplotlib.pyplot as plt
import random

#Global Varibales 
t = 1200 #mins
dt = 1 #time step, time in between frames
cyto_visc = 1.0 #normalized we can change
T = 60 #constant period, 60 minutes, typicall of cell-free Xenopus extract 
pressure_grad = 1.0 #normalized for now, will change later if neccesary 
vel_baseline = 1.0 #normalized but will change later
k = 1.0 #normalized but can change later
baseline_diffusion = 1.0
proportionality_const = 1.0


#Functions
def growth_rate(delta_rad):
    
    if delta_rad == 0:
        growth_rate = 0
    else:
        growth_rate = max(delta_rad,0) / dt
    
    return growth_rate

def highway_openness(rad,maximumradii):
    
    occ_fraction = (rad / maximumradii)
    
    if occ_fraction > 1:
        occ_fraction = 1
    return occ_fraction

def dyn_visc(occupationalfraction):
    
    dynamicvisco = cyto_visc * np.exp(-(k*occupationalfraction))
    
    return dynamicvisco

def local_flow_calc(vis, ugrowth, background):

    ulocal = background * vel_baseline
    ulocal += proportionality_const * (pressure_grad/vis)
    ulocal = max(ulocal, ugrowth) #boundary condition here:
        #it should never moved slower than the baseline OR its own boundary: ugrowth

    return ulocal

def eff_dif(visco):
    
    Deff = baseline_diffusion / visco
    
    if Deff < 1:
        Deff = 1
    
    return Deff

def transport_coe_calc(localvelocity, effectivediffusion):
    
    TransCoe = (localvelocity/vel_baseline) * (effectivediffusion/baseline_diffusion)
    
    return TransCoe

def delphaselagcalc(transportC):
    
    deltaphaselag =  (1 / transportC) * np.pi 
    #I think the max lag would be pi, this would be an inverse relationship. 
    #I doubt that the coupled relationship woul dbe able to sustain- 
    #if it started to lag behind more than 1/2 a wavelength. 
    
    return deltaphaselag


#arrays to store my data 
minute = []
curr_diameter = []
curr_radius = []
delta_radii = []
growth_r8 = []
highway_openness_fraction = []
dynamic_viscosity = []
gradientDV = []
local_vel = []
Effective_diffusion =[]
transport_coeffeient = []
delta_phase_lag = []


#delcaring to be safe! :) 
current_cycle = 0
cycle_time = 0


#Polyrate data, measured in ImageJ
#loop updating diameter based on the poly rates
for i in range(0, t, dt):
    
    background_transport = 1.0 + 0.15*np.sin(2*np.pi*i/(6*T)) 
    #accounts for background processes, like ATP, denisty,cyto org and the state
    
    current_cycle = i // T
    cycle_time = i % T
    max_rad = 1
    
    
    if current_cycle == 0:
        effective_diameter = 0.008
        effective_radius = effective_diameter / 2
        polyrate = 0
        max_rad = 1000 #no max bc no compartments
    
    else:
        
        #cycles 1-3   
        if 1 <= current_cycle < 3:
            max_rad = 507 + 20 #microns
            if cycle_time == 0:
                effective_diameter = 0.008
                effective_radius = effective_diameter/2
                polyrate = 0
            elif 1 <= cycle_time <= 34:
                polyrate = (10.583 * cycle_time - 68.421) * 2.6
                effective_diameter = polyrate
                if effective_diameter < 0.008:
                    effective_diameter = 0.008
                    effective_radius = 0.004
                else:
                    effective_diameter = polyrate #* random.uniform(.99, 1.01)
                    effective_radius = effective_diameter/2
                    
                    #restructuring, we see flux but the boundary is not clear...
                    #upon another look, I don't think it is 'shrinking' per say
                    #the highway areas do become crowded so I will add in a slight 
                    #growth factor 
            elif 34 < cycle_time < 60:
                polyrate = 1.015 * curr_diameter[-1] 
                effective_diameter = polyrate #* random.uniform(.99, 1.01)
                if effective_diameter < 0.008:
                    effective_diameter = 0.008
                    effective_radius = 0.004
                else:
                    effective_diameter = polyrate #* random.uniform(.99, 1.01)
                    effective_radius = effective_diameter/2
    
        #cycles 3-8
        elif 3 <= current_cycle < 8:
            max_rad = 266.5 + 20 #microns
            if cycle_time == 0:
                effective_diameter = 0.008
                effective_radius = effective_diameter/2
                polyrate = 0
            elif 1 <= cycle_time <= 42:
                polyrate = (7.6703 * cycle_time - 22.956) * 2.6
                if polyrate > 0:
                    if effective_diameter < 0.008:
                        effective_diameter = 0.008
                        effective_radius = 0.004
                    else:
                        effective_diameter = polyrate #* random.uniform(.99, 1.01)
                        effective_radius = effective_diameter/2
            elif 42 < cycle_time < 60:
                polyrate = (0.8247 * cycle_time + 174.46) * 2.6
                if effective_diameter < 0.008:
                    effective_diameter = 0.008
                    effective_radius = 0.004
                elif effective_diameter < polyrate: #shouldn't shrink af
                    effective_diameter = polyrate #* random.uniform(.99, 1.01)
                    effective_radius = effective_diameter/2
    
        #cyles 9 and past
        else:
            max_rad = 299 + 20 #microns
            if cycle_time == 0:
                effective_diameter = 0.008
                effective_radius = effective_diameter/2
                polyrate = 0
            elif 1 <= cycle_time <= 38:
                polyrate = ((6.7632 * cycle_time) - 6.5) * 2.6
                if effective_diameter < 0.008:
                    effective_diameter = 0.008
                    effective_radius = 0.004
                else:
                    effective_diameter = polyrate #* random.uniform(.99, 1.01)
                    effective_radius = effective_diameter/2
            elif 38 < cycle_time < 60:
                polyrate = (0.0296 * cycle_time + 212.39) * 2.6
                if effective_diameter < 0.008:
                    effective_diameter = 0.008
                elif effective_diameter < polyrate:
                    effective_diameter = polyrate #* random.uniform(.99, 1.01)
                    effective_radius = effective_diameter/2
    
    #adding data from if statements to arrays
    minute.append(i)
    
    if effective_radius > max_rad:
        effective_radius = max_rad
        effective_diameter = 2 * effective_radius
        
    curr_diameter.append(effective_diameter)
    curr_radius.append(effective_radius)

    #delta radii stuff
    if i == 0 or cycle_time == 0:
        delta_radii_carrier = 0
    else:
        delta_radii_carrier = curr_radius[-1] - curr_radius[-2]
    delta_radii.append(delta_radii_carrier)
    
    #growth rate stuff
    growth_rate_carrier = growth_rate(delta_radii_carrier)
    growth_r8.append(growth_rate_carrier)
    
    #highway fraction stuff
    highway_openness_carrier = highway_openness(effective_radius, max_rad)
    highway_openness_fraction.append(highway_openness_carrier)
    
    #dynamic viscoisty stuff
    dynamic_visocity_carrier = dyn_visc(highway_openness_carrier)
    dynamic_viscosity.append(dynamic_visocity_carrier)   
    
    #local Flow stuff
    ulocal_carrier = local_flow_calc(dynamic_visocity_carrier, growth_rate_carrier,background_transport)
    local_vel.append(ulocal_carrier)
    
    #effective diffusion stuff
    DEFF_carrier = eff_dif(dynamic_visocity_carrier)
    Effective_diffusion.append(DEFF_carrier)
    
    #Transport Coe stuff 
    transport_carrier = transport_coe_calc(ulocal_carrier, DEFF_carrier)
    transport_coeffeient.append(transport_carrier)
    
    #Change in Phase Lag 
    phase_lag_carrier = delphaselagcalc(transport_carrier)
    delta_phase_lag.append(phase_lag_carrier)
    
#END OF LOOP


#Printing stuff

for j in range(60,125,dt):
    print("MINUTE: ", minute[j], 
          "Radii: ", curr_radius[j],
          "Diameter: ", curr_diameter[j],
          "Change in radii: ", delta_radii[j],
          "Cluster Growth Velocity: ", growth_r8[j],
          "Highway Openness: ", highway_openness_fraction[j],
          "Dynamic Viscosity: ", dynamic_viscosity[j],
          "Local velocity: ", local_vel[j],
          "Effective Diffusion: ", Effective_diffusion[j],
          "Transport Coeffiecent: ",  transport_coeffeient[j],
          "Change in Phase Lag/Lead: ", delta_phase_lag[j])
    



plt.plot(minute, curr_diameter , label = 'diameter')
plt.legend() 
for x in range(60,1201,60):
    plt.axvline(x, linestyle='--', alpha=0.3, color='black') 
plt.show()

plt.plot(minute, curr_radius, label = 'radii')
plt.legend()  
for x in range(60,1201,60):
    plt.axvline(x, linestyle='--', alpha=0.3, color='black') 
plt.show()

plt.plot(minute, delta_radii, color = 'green', label = 'delta radii') 
plt.legend()  
for x in range(60,1201,60):
    plt.axvline(x, linestyle='--', alpha=0.3, color='black') 
plt.show()

plt.plot(minute, growth_r8, color = 'm', label = 'cluster growth velocity')
plt.legend()  
for x in range(60,1201,60):
    plt.axvline(x, linestyle='--', alpha=0.3, color='black') 
plt.show()

plt.plot(minute, highway_openness_fraction, color = 'blue', label = 'occupancy fraction')                                
plt.legend()  
for x in range(60,1201,60):
    plt.axvline(x, linestyle='--', alpha=0.3, color='black')                             
plt.show()

plt.plot(minute, dynamic_viscosity, color = 'darkorange', label = 'viscosity' )
plt.legend() 
for x in range(60,1201,60):
    plt.axvline(x, linestyle='--', alpha=0.3, color='black')  
plt.show()

plt.plot(minute, local_vel, color = 'red', label = 'local vel relative to baseline')
plt.legend()  
for x in range(60,1201,60):
    plt.axvline(x, linestyle='--', alpha=0.3, color='black') 
plt.show()

plt.plot(minute, Effective_diffusion, color = 'cyan', label = 'effective diffusion')
plt.legend()  
for x in range(60,1201,60):
    plt.axvline(x, linestyle='--', alpha=0.3, color='black') 
plt.show()
    
plt.plot(minute, transport_coeffeient, color = 'purple', label = 'transport coeffeient')
plt.legend()  
for x in range(60,1201,60):
    plt.axvline(x, linestyle='--', alpha=0.3, color='black') 
plt.show()

plt.plot(minute, delta_phase_lag, color = 'gray', label = 'change in phase lag')
plt.legend() 
y_ticks=(0, np.pi/2, np.pi)
#I think the max lag is about pi, not 2pi --> this assumption can be altered in later versions. 
y_labels = ['0','π/2','π']
plt.yticks(y_ticks, y_labels)
plt.xlabel("Time (minutes)")
plt.ylabel("Relative Phase Lag (radians)")
for x in range(60,1201,60):
    plt.axvline(x, linestyle='--', alpha=0.3, color='black')  
plt.show()

     