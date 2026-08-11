#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 30 17:20:58 2026

@author: maggie
"""

#global libraries
import numpy as np
import matplotlib.pyplot as plt
import random
from scipy import stats

#global variables 
t = 1200 #mins
dt = 1 #timestep
k_eff = 1.0e-6 #kind of like our spring constant 
rest_radii = 10 #can change later, but essentially our 'springs' rest radii 
gamma = 1.0 #scaling constant for energy to force
max_force_releasable = 0.3 #Newtons, adjust this to match the largest force produced by model sim
baseline_CDK1_wave_speed = 10 #microns / second
baseline_diffusion = 1.0 
drag_force = 0.05 #edit later depending on order of mag
tau = 1.0
alpha = 2175 #experimental data 
beta = 1
theta = 0.65 # will tweak this later, see living doc for why this is the number it is.... 
scaling_const = 3e3 # gets us onto the order or magitude same as the baseline
scaling_std = 1805 #see excel
lag_mins_carrier = 8 #starting lag, experimental data 


#Functions
def growth_rate(delta_rad):
    
    if delta_rad == 0:
        growth_rate = 0
    else:
        growth_rate = max(delta_rad,0) / dt
    
    return growth_rate

def last_nonzero(arr):
    for i in range(len(arr)-1, -1, -1):
        if arr[i] > 0:
            print(f"Found energy {arr[i]:.6f} at index {i}")
            return arr[i]
    return 0

def calculated_mechanical_energy(radii_parameter):

    stored_mechanical_energy = (0.5 * k_eff * (radii_parameter - rest_radii)**2) 
    
    return stored_mechanical_energy

def calculated_transport_coe(wavespeed,diffusion):
    
    Transport_COE = (wavespeed/baseline_CDK1_wave_speed) * (diffusion/baseline_diffusion)
    
    return Transport_COE

def calculated_phase_shift(transportcoe):
    
    deltaphaseshift = (1 / transportcoe)  # added to make look better but should find reason as to why... 
    
    return  deltaphaseshift



#initalizing 
cycle_time = 0
current_cycle = 1
Period = 46
mixing_val = 0.5 #consider this our baseline, mixing val is between 0 and 1. 
initial_CDK1_wave_speed = baseline_CDK1_wave_speed
initial_diffusion = baseline_diffusion
force_released = 0
prev_force_release = 0.04 #assuming this largest force, itital compartment formation and depoly
current_stored_energy = 0

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
initial_CDK1_wavespeed_arr = []
initial_diffusion_arr = []
curr_wave_speed = []
curr_diffusion_rate = []
Transport_coe_arr = []
delta_phase_shift = []
cycle_periods = []
cycle_end_times = []
initial_mixing_value = []
released_force = []
inital_wavespeed = []
inital_dif = []
mixing_val_arr_test = []
lag_mins = []
cycle_lags = []
cycle_numbers = []

released_force.append(prev_force_release)


for i in range(0, t, dt):
     
    cycle_time += dt #cycle_time = time since last wave
    
    
    if current_cycle == 0:
        if cycle_time > 0:
            growth_slope = 6.5088
            grow_time = Period
            polyrate = (growth_slope * cycle_time) + 178.63
            if polyrate > 0:
                effective_diameter = polyrate * 2.6
                effective_radius = effective_diameter / 2
            else:
                 effective_diameter = 0.008
                 effective_radius = 0.004
            max_diameter = ((1.1517 *( growth_slope**2)) - (11.414 * growth_slope) + 586.72) * 2.6
            max_radii = max_diameter / 2
    
    elif current_cycle == 1:
        if cycle_time > 0:
            growth_slope = 10.005
            grow_time = Period
            polyrate = (growth_slope * cycle_time) - 9.5344
            if polyrate > 0:
                effective_diameter = polyrate * 2.6
                effective_radius = effective_diameter / 2
            else:
                 effective_diameter = 0.008
                 effective_radius = 0.004
            max_diameter = ((1.1517 *( growth_slope**2)) - (11.414 * growth_slope) + 586.72) * 2.6
            max_radii = max_diameter / 2
    
    elif current_cycle == 2:
        if 0 < cycle_time <= 20:
            growth_slope = 15.155
            grow_time = 20
            polyrate = (growth_slope * cycle_time) - 33.545
            if polyrate > 0:
                effective_diameter = polyrate * 2.6
                effective_radius = effective_diameter / 2
            else:
                 effective_diameter = 0.008
                 effective_radius = 0.004
            max_diameter = ((1.1517 *( growth_slope**2)) - (11.414 * growth_slope) + 586.72) * 2.6
            max_radii = max_diameter / 2
        else:
            stablization_slope = 2.1033
            polyrate = (stablization_slope * cycle_time) + 212.53
            effective_diameter = polyrate * 2.6
            effective_radius = effective_diameter / 2 
   
    elif current_cycle == 3:
        if 0 < cycle_time <= 12:
            growth_slope = 23.214
            grow_time = 12
            polyrate = (growth_slope * cycle_time) - 47.143
            if polyrate > 0:
                effective_diameter = polyrate * 2.6
                effective_radius = effective_diameter / 2
            else:
                 effective_diameter = 0.008
                 effective_radius = 0.004
            max_diameter = ((1.1517 *( growth_slope**2)) - (11.414 * growth_slope) + 586.72) * 2.6
            max_radii = max_diameter / 2
        else:
            stablization_slope = 1.2946
            polyrate = (stablization_slope * cycle_time) + 282.75
            effective_diameter = polyrate * 2.6
            effective_radius = effective_diameter / 2
    
    elif current_cycle == 4:
        if 0 < cycle_time <= 16:
            growth_slope = 22.55
            grow_time = 16
            polyrate = (growth_slope * cycle_time) + 11.933
            if polyrate > 0:
                effective_diameter = polyrate * 2.6
                effective_radius = effective_diameter / 2
            else:
                 effective_diameter = 0.008
                 effective_radius = 0.004
            max_diameter = ((1.1517 *( growth_slope**2)) - (11.414 * growth_slope) + 586.72) * 2.6
            max_radii = max_diameter / 2
        else:
            stablization_slope = -0.8636
            polyrate = (stablization_slope * cycle_time) + 369.18
            effective_diameter = polyrate * 2.6
            effective_radius = effective_diameter / 2
    
    elif current_cycle == 5:
        if 0 < cycle_time <= 22:
            growth_slope = 13.357
            grow_time = 22
            polyrate = (growth_slope * cycle_time) - 46.923
            if polyrate > 0:
                effective_diameter = polyrate * 2.6
                effective_radius = effective_diameter / 2
            else:
                effective_diameter = 0.008
                effective_radius = 0.004
            max_diameter = ((1.1517 *( growth_slope**2)) - (11.414 * growth_slope) + 586.72) * 2.6
            max_radii = max_diameter / 2
        else:
            stablization_slope = -0.0824
            polyrate = (stablization_slope * cycle_time) + 278.74
            effective_diameter = polyrate * 2.6
            effective_radius = effective_diameter / 2
    
    elif current_cycle == 6:
        if 0 < cycle_time <= 16:
            growth_slope = 13.083
            grow_time = 16
            polyrate = (growth_slope * cycle_time) - 16.889
            if polyrate > 0:
                effective_diameter = polyrate * 2.6
                effective_radius = effective_diameter / 2
            else:
                 effective_diameter = 0.008
                 effective_radius = 0.004
            max_diameter = ((1.1517 *( growth_slope**2)) - (11.414 * growth_slope) + 586.72) * 2.6
            max_radii = max_diameter / 2
        else:
            stablization_slope = -1.4118
            polyrate = (stablization_slope * cycle_time) + 217.21
            effective_diameter = polyrate * 2.6
            effective_radius = effective_diameter / 2
        
    elif current_cycle == 7:
        if 0 < cycle_time <= 14:
            growth_slope = 13.512
            grow_time = 14
            polyrate = (growth_slope * cycle_time) + 2.9167
            if polyrate > 0:
                effective_diameter = polyrate * 2.6
                effective_radius = effective_diameter / 2
            else:
                 effective_diameter = 0.008
                 effective_radius = 0.004
            max_diameter = ((1.1517 *( growth_slope**2)) - (11.414 * growth_slope) + 586.72) * 2.6
            max_radii = max_diameter / 2
        else:
          stablization_slope = -0.2632
          polyrate = (stablization_slope * cycle_time) + 194.24
          effective_diameter = polyrate * 2.6
          effective_radius = effective_diameter / 2  
    
    elif current_cycle == 8:
        if 0 < cycle_time <= 20:
            growth_slope =  9.6
            grow_time = 20
            polyrate = (growth_slope * cycle_time) - 14.727
            if polyrate > 0:
                effective_diameter = polyrate * 2.6
                effective_radius = effective_diameter / 2
            else:
                 effective_diameter = 0.008
                 effective_radius = 0.004
            max_diameter = ((1.1517 *( growth_slope**2)) - (11.414 * growth_slope) + 586.72) * 2.6
            max_radii = max_diameter / 2
        else:
          stablization_slope = 0.3137
          polyrate = (stablization_slope * cycle_time) + 155.33
          effective_diameter = polyrate * 2.6
          effective_radius = effective_diameter / 2 
    
    elif current_cycle == 9:
        if 0 < cycle_time <= 28:
            growth_slope = 7.1071
            grow_time = 28
            polyrate = (growth_slope * cycle_time) - 33.833
            if polyrate > 0:
               effective_diameter = polyrate * 2.6
               effective_radius = effective_diameter / 2
            else:
               effective_diameter = 0.008
               effective_radius = 0.004
            max_diameter = ((1.1517 *( growth_slope**2)) - (11.414 * growth_slope) + 586.72) * 2.6
            max_radii = max_diameter / 2
        else:
           stablization_slope =  0.3827
           polyrate = (stablization_slope * cycle_time) + 146.73
           effective_diameter = polyrate * 2.6
           effective_radius = effective_diameter / 2  
    
    elif current_cycle == 10:
        if 0 < cycle_time <= 28:
            growth_slope = 6.3482
            grow_time = 28
            polyrate = (growth_slope * cycle_time) - 16.875
            if polyrate > 0:
               effective_diameter = polyrate * 2.6
               effective_radius = effective_diameter / 2
            else:
                effective_diameter = 0.008
                effective_radius = 0.004
            max_diameter = ((1.1517 *( growth_slope**2)) - (11.414 * growth_slope) + 586.72) * 2.6
            max_radii = max_diameter / 2
        else:
            stablization_slope =  0.3111
            polyrate = (stablization_slope * cycle_time) + 162.5
            effective_diameter = polyrate * 2.6
            effective_radius = effective_diameter / 2  
    
    elif current_cycle == 11:
        if 0 < cycle_time <= 34:
            growth_slope = 5.485
            grow_time = 34
            polyrate = (growth_slope * cycle_time) - 22.69
            if polyrate > 0:
               effective_diameter = polyrate * 2.6
               effective_radius = effective_diameter / 2
            else:
                effective_diameter = 0.008
                effective_radius = 0.004
            max_diameter = ((1.1517 *( growth_slope**2)) - (11.414 * growth_slope) + 586.72) * 2.6
            max_radii = max_diameter / 2
        else:
            stablization_slope =  0.1921
            polyrate = (stablization_slope * cycle_time) + 146.67
            effective_diameter = polyrate * 2.6
            effective_radius = effective_diameter / 2  
    
    elif current_cycle == 12:
        if 0 < cycle_time <= 36: 
            growth_slope = 5.4123
            grow_time = 36
            polyrate = (growth_slope * cycle_time) - 16.368
            if polyrate > 0:
               effective_diameter = polyrate * 2.6
               effective_radius = effective_diameter / 2
            else:
                effective_diameter = 0.008
                effective_radius = 0.004
            max_diameter = ((1.1517 *( growth_slope**2)) - (11.414 * growth_slope) + 586.72) * 2.6
            max_radii = max_diameter / 2
        else:
           stablization_slope =  0.1597
           polyrate = (stablization_slope * cycle_time) + 149.22
           effective_diameter = polyrate * 2.6
           effective_radius = effective_diameter / 2 
         
    
    #adding data from if statements to arrays
    minute.append(i)
      
    if effective_radius > max_radii:
        effective_radius = max_radii
        effective_diameter = 2 * effective_radius
         
    curr_diameter.append(effective_diameter)
    curr_radius.append(effective_radius)

    #delta radii stuff
    if i == 0 or cycle_time == 1:
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
    if growth_rate_carrier > 0:
        current_stored_energy = calculated_mechanical_energy(effective_radius)

    stored_mechanical_energy.append(current_stored_energy)
    #have to work in a catch here to catch when mechanical energy = 0 due to
    #the growth rate being 0, mechanical energy should remain stored regardless

    
    #CDK1 wave hit loop check 
    if cycle_time >= Period:
        print("CDK1 wave hit")
        
        cycle_periods.append(cycle_time)
        cycle_end_times.append(minute[i]+1)
        
        cycle_lags.append(lag_mins_carrier)
        cycle_numbers.append(current_cycle)
        
        #force calc when the wave hits the cluster
        force_released = gamma * current_stored_energy
        released_force.append(force_released) #issue, this is changing when mixing does
        #the force released is independent of mixing. 
        
        mixing_val = prev_force_release / max_force_releasable
        
        #we are looking at how the pervious force affected net mixing
        prev_force_release = force_released
    
        initial_mixing_value.append(mixing_val)
                                                            
        mixing_val = max(mixing_val, 1e-6)
        
        
        #need to work on this scaling, might make alpha variable
        #should not be so many periods...
        #using the excel relation between alpha and cycle number 
        #maybe alter mixing coe 
        #also maybe some sort of time constraint, cycle number, there is an exponential trend 
        #end should exponentially trend upwards, 
        
        Period = 46 + (alpha  * ( growth_slope**-0.4206) * (max_diameter**-0.56135) * (mixing_val**-theta))
        # 0.05 is an assumption, makes the plot look good lol
        
        initial_CDK1_wave_speed = ((.532 / Period) * scaling_const) * (mixing_val+1) 
        #532 is a wavelength approx
        #wavespeed is affected by force, whihc affects mixing
        inital_wavespeed.append(initial_CDK1_wave_speed)
        
        
        initial_diffusion = (mixing_val*baseline_diffusion) + baseline_diffusion
        inital_dif.append(initial_diffusion)
        
        current_cycle += 1
        cycle_time = 0
        
        effective_radius = 0.004
        effective_diameter = 0.008
        
        
    
    #may not need to do this right now, itial wavespeed might be enough
    wave_speed_decay = max(baseline_CDK1_wave_speed , (initial_CDK1_wave_speed) - (drag_force * cycle_time))
    curr_wave_speed.append(wave_speed_decay)
     
    diffusion_deacy = max(baseline_diffusion,(baseline_diffusion + (initial_diffusion-baseline_diffusion)*np.exp(-cycle_time/Period)))
    curr_diffusion_rate.append(diffusion_deacy)
     
    Transport_coe = calculated_transport_coe(wave_speed_decay, diffusion_deacy)
    Transport_coe_arr.append(Transport_coe)
     
    delta_phase_shift_carrier = calculated_phase_shift(Transport_coe)
    delta_phase_shift.append(delta_phase_shift_carrier)
    
    lag_mins_carrier = Period * (8/46) * delta_phase_shift_carrier
    #(8/46) is the inital lag period. 8 minutes for a 46 min cycle.
    lag_mins.append(lag_mins_carrier)
    

    T.append(Period)
    mixing_val_arr.append(mixing_val)
    force_released_arr.append(force_released)
    initial_CDK1_wavespeed_arr.append(initial_CDK1_wave_speed)
    initial_diffusion_arr.append(initial_diffusion)

 #END OF LOOP ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

released_force.pop()

  #printing stuff
for j in range(45,47,dt):
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
          "Initial Diffusion: ", initial_diffusion_arr[j])
     
    
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

plt.plot(cycle_end_times, released_force , marker='o', label = 'Force Released', color = 'darkorange')
plt.legend()  
plt.show()

plt.plot(cycle_end_times, initial_mixing_value , marker='o', label = 'Mixing Induced by Force', color = 'pink')
plt.legend()  
plt.show()

plt.plot(cycle_end_times, cycle_periods, marker='o', color='black', label='Cycle Period')
plt.xlabel("Time (min)")
plt.ylabel("Period (min)")
plt.axvline(952, linestyle='--', alpha=0.3, color='black') 
plt.legend()
plt.show()

plt.plot(cycle_end_times, inital_wavespeed , marker='o', label = 'Predicted Initial Wavespeed after Force Remodeling', color = 'green')
plt.legend()  
plt.show()


plt.plot(minute, curr_wave_speed , label = 'Predicted Wavespeed over time', color = 'blue')
plt.xlabel("Time (min)")
plt.ylabel("Wave Speed (microns/min)")
plt.legend()  
plt.show()       


plt.plot(cycle_end_times, inital_dif , marker='o', label = 'Predicted Initial Diffusion after Force Remodeling', color = 'red')
plt.legend()  
plt.show()


plt.plot(minute, curr_diffusion_rate , label = 'Predicted rate of Diffusion over Time')
plt.legend()  
plt.show() 

plt.plot(minute, Transport_coe_arr , label = 'Transport Coeffiecent over Time', color = 'purple')
plt.legend()  
plt.show()   


plt.plot(minute, delta_phase_shift, label = 'Phase Shift over Time', color = 'green')
y_ticks=(0, np.pi, 2 * np.pi)
y_labels = ['0','π','2π']
plt.yticks(y_ticks, y_labels)
plt.axvline(952, linestyle='--', alpha=0.3, color='black') 
plt.axvline(46, linestyle='--', alpha=0.3, color='black') 
plt.legend()  
plt.show()  

plt.plot(minute, lag_mins, label = 'minute lag overtime', color = 'm')  
plt.legend()  
plt.show() 


import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable

x = np.arange(len(cycle_periods))

slope, intercept, r_value, p_value, std_err = stats.linregress(x, cycle_periods)

# 5. Create and style the plot
plt.figure(figsize=(9, 7))

plt.scatter(x, cycle_periods, color="blue", s=80, zorder=3, vmin=-2)
plt.plot(x, slope*x + intercept, color="black", linewidth=2, label='Linear Fit')

equation = f"y = {slope:.2f}x + {intercept:.2f}\n$R^2$ = {r_value**2:.2f}"

plt.text(
    0.05, 0.95,
    equation,
    transform=plt.gca().transAxes,
    fontsize=10,
    verticalalignment='top',
    bbox=dict(facecolor='white', alpha=0.8)
)

plt.ylim(40, 140)
plt.title('Period Length vs Cycle Number (Model Results)', fontsize=22, fontweight='bold', pad=15)
plt.xlabel('Cycle #', fontsize=16, labelpad=8)
plt.ylabel('Period Length (mins)', fontsize=16, labelpad=8)
plt.grid(True, linestyle='--', alpha=0.6)
plt.tight_layout()
plt.show()



plt.figure(figsize=(9, 7))
slope, intercept, r_value, p_value, std_err = stats.linregress(x, cycle_lags)
plt.plot(x, slope*x + intercept, color="black")
equation = f"y = {slope:.2f}x + {intercept:.2f}\n$R^2$ = {r_value**2:.2f}"
plt.scatter(x, cycle_lags, marker='o', color = "blue")

plt.text(
    0.05, 0.95,
    equation,
    transform=plt.gca().transAxes,
    fontsize=10,
    verticalalignment='top',
    bbox=dict(facecolor='white', alpha=0.8)
)

plt.xlabel("Cycle #", fontsize=16, labelpad=8)
plt.ylabel("Phase Lag (min)", fontsize=16, labelpad=8)
plt.title("Lag Between Activity Per Cycle (Model Results)",fontsize=22, fontweight='bold', pad=15)
plt.grid(True, linestyle='--', alpha=0.6)
plt.tight_layout()
plt.show()







  

    