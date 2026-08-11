#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 29 13:36:06 2026

@author: maggie
"""

#global libraries
import numpy as np
import matplotlib.pyplot as plt

#global variables 
t = 1200 #mins
dt = 1 #timestep
k_eff = 1.0e-6 #kind of like our spring constant 
rest_radii = 10 #can change later, but essentially our 'springs' rest radii 
gamma = 1.0 #scaling constant for energy to force
max_force_releasable = 0.14 #Newtons, adjust this to match the largest force later
baseline_crowding = 0.5 #start in the middle 
baseline_CDK1_wave_speed = 1.0 #microns / second
baseline_diffusion = 1.0 
drag_force = 0.05 #edit later depending on order of mag
tau = 1.0


#Functions
def growth_rate(delta_rad):
    
    if delta_rad == 0:
        growth_rate = 0
    else:
        growth_rate = max(delta_rad,0) / dt
    
    return growth_rate

def calculated_mechanical_energy(radii_parameter):
    
    stored_mechanical_energy = 0.5 * k_eff * (radii_parameter - rest_radii)**2
    
    return stored_mechanical_energy

def calculated_transport_coe(wavespeed,diffusion):
    
    Transport_COE = (wavespeed/baseline_CDK1_wave_speed) * (diffusion/baseline_diffusion)
    
    return Transport_COE

def calculated_phase_shift(transportcoe):
    
    deltaphaseshift = 1 / transportcoe
    
    return  deltaphaseshift

#storage arrays
minute = []
curr_radius = []
curr_diameter = []
delta_radii = []
growth_r8 = []
stored_mechanical_energy = []
force_released_arr = []
mixing_val_arr = []
T = []
crowding_coe_arr = []
initial_CDK1_wavespeed_arr = []
initial_diffusion_arr = []
curr_wave_speed = []
curr_diffusion_rate = []
Transport_coe_arr = []
delta_phase_shift = []

#initalizing 
cycle_time = 0
current_cycle = 1
Period = 75 # starting period length, minutes, from current data 
mixing_val = 0.5 #consider this our baseline
initial_CDK1_wave_speed = baseline_CDK1_wave_speed
initial_diffusion = baseline_diffusion
crowding_coe = 1.0
force_released = 0

for i in range(0, t, dt):
     
    cycle_time += dt #cycle_time = time since last wave
    
    
    if current_cycle == 0:
        effective_diameter = 0.008
        effective_radius = effective_diameter / 2
        polyrate = 0
        max_rad = 1e10 #no max bc no compartments
    
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
            elif 34 < cycle_time:
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
            max_rad = 273 + 100 #microns
            if cycle_time == 0:
                effective_diameter = 0.008
                effective_radius = effective_diameter/2
                polyrate = 0
            elif 1 <= cycle_time <= 42:
                polyrate = (7.6703 * cycle_time - 22.956) * 2.6
                #effective_diameter = polyrate
                if polyrate > 0:
                    if effective_diameter < 0.008:
                        effective_diameter = 0.008
                        effective_radius = 0.004
                    else:
                        effective_diameter = polyrate #* random.uniform(.99, 1.01)
                        effective_radius = effective_diameter/2
            elif 42 < cycle_time:
                polyrate = (0.8247 * cycle_time + 174.46) * 2.6
                effective_diameter = polyrate
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
                effective_diameter = polyrate
                if effective_diameter < 0.008:
                    effective_diameter = 0.008
                    effective_radius = 0.004
                else:
                    effective_diameter = polyrate #* random.uniform(.99, 1.01)
                    effective_radius = effective_diameter/2
            elif 38 < cycle_time:
                polyrate = (0.0296 * cycle_time + 212.39) * 2.6
                effective_diameter = polyrate
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
         delta_radii_carrier = curr_radius[i] - curr_radius[i-1]
         if delta_radii_carrier < 0:
             delta_radii_carrier = 0
     
    delta_radii.append(delta_radii_carrier)
        
    #growth rate stuff
    growth_rate_carrier = growth_rate(delta_radii_carrier)
    growth_r8.append(growth_rate_carrier)
    
    #stored mechanical energy stuff 
    mechanical_energy_carrier = calculated_mechanical_energy(effective_radius)
    stored_mechanical_energy.append(mechanical_energy_carrier)
    
    if cycle_time < Period:
        force_released = 0
    
    #CDK1 wave hit loop check 
    if cycle_time >= Period:
        print("CDK1 wave hit")
        
        #force calc when the wave hits the cluster
        force_released = gamma * stored_mechanical_energy[i]
        
        mixing_val = force_released / max_force_releasable #want to normalize this to 1,
                                                            #hence prev line
        mixing_val = max(mixing_val, 1e-6)
        
        if mixing_val > 0.5:
            Period = 75 - 12*(1-mixing_val)
        else:
            Period = 75 + 12*(1-mixing_val)
        
        crowding_coe = baseline_crowding / (mixing_val)
        
        initial_CDK1_wave_speed = baseline_CDK1_wave_speed / crowding_coe
        
        initial_diffusion = (mixing_val*baseline_diffusion) + baseline_diffusion
        
        current_cycle += 1
        cycle_time = 0
        
        effective_radius = 0.004
        effective_diameter = 0.008
    
    wave_speed_decay = max(baseline_CDK1_wave_speed, (initial_CDK1_wave_speed) - (drag_force * cycle_time))#this is gonna throw an error
    curr_wave_speed.append(wave_speed_decay)
    
    diffusion_deacy = max(baseline_diffusion,(baseline_diffusion + (initial_diffusion-baseline_diffusion)*np.exp(-cycle_time/tau)))#i think this is also gonna throw an error
    curr_diffusion_rate.append(diffusion_deacy)
    
    Transport_coe = calculated_transport_coe(wave_speed_decay, diffusion_deacy)
    Transport_coe_arr.append(Transport_coe)
    
    delta_phase_shift_carrier = calculated_phase_shift(Transport_coe)
    delta_phase_shift.append(delta_phase_shift_carrier)


    T.append(Period)
    crowding_coe_arr.append((crowding_coe))
    mixing_val_arr.append(mixing_val)
    force_released_arr.append(force_released)
    initial_CDK1_wavespeed_arr.append(initial_CDK1_wave_speed)
    initial_diffusion_arr.append(initial_diffusion)

#END OF LOOP ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

 #printing stuff
for j in range(0,t,20):
     print("MINUTE: ", minute[j], 
           "Radii: ", curr_radius[j],
           "Diameter: ", curr_diameter[j],
           "Change in radii: ", delta_radii[j],
           "Rate of Growth: ", growth_r8[j], #this currently unused....
           "Stored Mechanical Energy: ", stored_mechanical_energy[j],
           "Current Wave Speed: ", curr_wave_speed[j],
           "Current Rate of Diffusion: ", curr_diffusion_rate[j],
           "Transport Coefficent: ", Transport_coe_arr[j],
           "Change in Phase Shift: ", delta_phase_shift[j],
           "Period: ", T[j],
           "Mixing Coeifficent: ", mixing_val_arr[j],
           "Initial Wavespeed: ", initial_CDK1_wavespeed_arr[j],  
          "Initial Diffusion: ", initial_diffusion_arr[j],
          "Crowding Coefficent: ", crowding_coe_arr[j])
    
        
#plots
plt.plot(minute, curr_diameter , label = 'diameter')
plt.legend() 
plt.show() #this is a little wonky but radii seems fine

plt.plot(minute, curr_radius, label = 'radii')
plt.legend()  
plt.show()  

plt.plot(minute, delta_radii , label = 'Delta Radii', color = 'red')
plt.legend()  
plt.show()       

plt.plot(minute, growth_r8 , label = 'Growth Rate', color = 'green')
plt.legend() #same as delta_radii 
plt.show()

plt.plot(minute, stored_mechanical_energy , label = 'Stored Mechanical Energy', color = 'blue')
plt.legend()  
plt.show()

plt.plot(minute, force_released_arr , label = 'Force Released', color = 'darkorange')
plt.legend()  
plt.show() 

plt.plot(minute, T, color='black', linewidth=2, label='Oscillation Period')
plt.legend()
plt.show()


 
