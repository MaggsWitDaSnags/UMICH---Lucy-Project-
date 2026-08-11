#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 17 11:51:27 2026

@author: maggie
"""
#POS 3, compartment
import numpy as np
import random

#global variable declarations
max_volume = (0.5 * np.pi * ((520)**3)) # (microns^3)
t = 1200 # mins
dt = 2 # mins 
cyto_visc = 1.20e-03 # g/micron*min
cyto_diff = 5 * 60 # microns^2/min

variance = random.randint(0, 15)
sign = random.choice([-1, 1])
T = 60

k = 7

#function definitions
#math is iffy but it is accurate in terms of what trends we see experimentally 
def volume(dia):
    radius = 0.5 * dia
    vol = (4/3) * np.pi * (radius**3)
    return vol

def local_density(curr_vol):
    if curr_vol >= max_volume:
        return 1.0
    else:
        local_den = curr_vol / max_volume
        return local_den

def dynamic_visc(localden):
    dyn_visc = np.exp(k * localden)
    return dyn_visc

def eff_dif(dynamicvisco):
    eff_diffusion = cyto_diff / dynamicvisco
    return eff_diffusion 

def wavespeed(diffusion):
    c = 2 * np.sqrt(diffusion)
    return c


history = {
    "minute": [],
    "diameter": [],
    "volume": [],
    "density": [],
    "viscosity": [],
    "diffusion": [],
    "wavespeed": [],
    "avg_fluid_vel": []  # New metric tracked dynamically
}

effective_diameter = 0.001 #microns, simple tubulin

current_cycle = 0
cycle_time = 0


#iterate through all 1200 minutes
for i in range(0, t, 1):
    current_cycle = i // T
    cycle_time = i % T
    
    #cycle 0, pre polyermization 
    if current_cycle == 0:
        effective_diameter = 0.001
        polyrate = 0
   
    else:
        
        #cycles 1-3   
        if 1 <= current_cycle < 3:
            if cycle_time == 0:
                effective_diameter = 0.001
                polyrate = 0
            elif 1 <= cycle_time <= 35:
                polyrate = (10.583 * cycle_time - 68.421) * 2.6
                if polyrate > 0:
                    effective_diameter = polyrate * random.uniform(.9, 1.1)
            elif 35 < cycle_time:
                polyrate = (-5.3734 * cycle_time + 352.23) * 2.6
                if polyrate > 0:
                    effective_diameter = polyrate * random.uniform(.9, 1.1)
    
        #cycles 3-8
        elif 3 <= current_cycle < 9:
            if cycle_time == 0:
                effective_diameter = 0.001
                polyrate = 0
            elif 1 <= cycle_time <= 24:
                polyrate = (7.6703 * cycle_time - 22.956) * 2.6
                if polyrate > 0:
                    effective_diameter = polyrate * random.uniform(.9, 1.1)
            elif 24 < cycle_time:
                polyrate = (0.8247 * cycle_time + 174.46) * 2.6
                if polyrate > 0:
                    effective_diameter = polyrate * random.uniform(.9, 1.1)
    
        #cyles 9 and past
        else:
            if cycle_time == 0:
                effective_diameter = 0.001
                polyrate = 0
            elif 1 <= cycle_time <= 38:
                polyrate = (6.7632 * cycle_time - 6.5) * 2.6
                if polyrate > 0:
                    effective_diameter = polyrate * random.uniform(.9, 1.1)
            elif 38 < cycle_time:
                polyrate = (0.0296 * cycle_time + 212.39) * 2.6
                if polyrate > 0:
                    effective_diameter = polyrate * random.uniform(.9, 1.1)
                
    #calculations
    current_volume = volume(effective_diameter)
    l_density = local_density(current_volume)
    #if period is too long, we start to see crosslinking which in turn affects the local denisty 
    if cycle_time > 60: 
        l_density = l_density * (cycle_time/60)
    dynamic_viscosity = dynamic_visc(l_density)
    effective_diffusion = eff_dif(dynamic_viscosity)
    wave_speed = wavespeed(effective_diffusion)
    
    #Logging the data to our history 
    history["minute"].append(i)
    history["diameter"].append(effective_diameter)
    history["volume"].append(current_volume)
    history["density"].append(l_density)
    history["viscosity"].append(dynamic_viscosity)
    history["diffusion"].append(effective_diffusion)
    history["wavespeed"].append(wave_speed)
    
    cycle_time += 1

total_cycles = t // T
cycle_starts = [c * T for c in range(total_cycles)]
cycle_ends = [(c + 1) * T - 1 for c in range(total_cycles)]

#~~~~~~~~NEW STUFF~~~~~~~~~~~~~~~
#THIS IS WAVESPEED STUFF SPECIFIC TO LOCATION
#wave location, lateral prop
length = 1250 
step = 1
tubulin_core_loc = 625 

#making our grid
x_grid = np.arange(0, length+ step, step)

history_nopoly = {
    "minutes": [],
    "wavespeed": []
}

history_poly = {
    "minutes": [],
    "wavespeed": []
}

current_cycle = 0
cycle_time = 0


for i in range(0, t, 1):
    
    #determining which cycle we are in dependning on how much time has elapsed
    current_cycle = i // T
    cycle_time = i % T
        
    base_diff = history["diffusion"][i]
    current_volume = history["volume"][i]

    #in the grid we want to measure how the wavespeed is affected in unploymerized area vs polyermized 

    #for the spaces where there are NO POLYMERIZATION
    highway_density = np.zeros_like(x_grid)
    highway_diffusion_coe = eff_dif(1.0) #empty cyto
    highway_wavespeed = wavespeed(highway_diffusion_coe)
    highway_transit_time = length/highway_wavespeed

    
    #POLYMERIZATION AREAS
    cluster_density = 0.8

    #we are making a guassian 'plot' of the denisty since it is not unifrom
    #most dense at center, one point
    cluster_radius = history["diameter"][i] / 2.0
    decay_width = max(cluster_radius, 10.0)
    cluster_density_profile = 0.8 * np.exp(-0.5 * ((x_grid - tubulin_core_loc) / decay_width)**2)
    
    #normal calcs 
    center_idx = len(x_grid) // 2
    cluster_visc = dynamic_visc(cluster_density_profile)
    cluster_dif_coe = eff_dif(cluster_visc)
    cluster_wavespeed = wavespeed(base_diff / cluster_visc)
    cluster_transit_time = np.sum(step / cluster_wavespeed)
    
    avg_highway_ws = np.mean(highway_wavespeed)
    avg_cluster_ws = np.mean(cluster_wavespeed)
    
    #finding wavespeed at center of graph 
    history_nopoly["minutes"].append(highway_transit_time)
    history_nopoly["wavespeed"].append(avg_highway_ws)
    history_poly["minutes"].append(cluster_transit_time)
    history_poly["wavespeed"].append(avg_cluster_ws)
    
    # --- Stokes Flow Dynamic Core Calculations ---
    peak_p = current_volume / max_volume
    
    if peak_p == 0 or history["diameter"][i] <= 0.001:
        history["avg_fluid_vel"].append(0.0)
    else:
        # Scale local pressure drop profile across space matching your Gaussian curve
        local_pressure_profile = cluster_density_profile * peak_p
        dp_dx_cluster = np.gradient(local_pressure_profile, step)
        
        # Determine velocity integration vectors balancing viscosity transformations
        d2u_dx2_cluster = dp_dx_cluster / cluster_visc
        u_cluster = np.cumsum(np.cumsum(d2u_dx2_cluster) * step) * step
        
        avg_fluid_velocity = np.mean(np.abs(u_cluster))
        # Zero catch safeguard for microscale computational noise limits
        if avg_fluid_velocity < 1e-12:
            history["avg_fluid_vel"].append(0.0)
        else:
            history["avg_fluid_vel"].append(avg_fluid_velocity)

    cycle_time += 1
#~~~~~~~ end of new stuff ~~~~~~~~~~~~~~


#printing the data
print(f"\n========================== DETAILED 10-MINUTE INTERVALS (DYNAMIC CYCLES) ==========================")
for c in range(total_cycles):
    cycle_start_min = cycle_starts[c]
    cycle_end_min = cycle_ends[c]
    if cycle_start_min >= len(history["minute"]):
        break
    if cycle_end_min >= len(history["minute"]):
        cycle_end_min = len(history["minute"]) - 1    
    cycle_length = cycle_end_min - cycle_start_min + 1
    
    print(f"\n--- CYCLE {c} (Minutes {cycle_start_min} to {cycle_end_min}, Length = {cycle_length} mins) ---")
    print(f"{'Time in Cycle':<15} | {'Global Min':<11} | {'Diameter (µm)':<15} | {'Density':<10} | {'Wavespeed (µm/min)':<18} | {'Avg Fluid Vel':<15}")
    print("-" * 101)
    
    # Store dynamic tracker to find highest entry in current loop execution block
    max_vel_in_cycle = 0.0
    
    for cycle_time in range(0, cycle_length, 10):
        global_min = cycle_start_min + cycle_time
        dia = history["diameter"][global_min]
        den = history["density"][global_min]
        ws = history["wavespeed"][global_min]
        v_fluid = history["avg_fluid_vel"][global_min]
        
        # Track highest velocity value
        if v_fluid > max_vel_in_cycle:
            max_vel_in_cycle = v_fluid
            
        # Format string check handles zero catch values natively 
        v_fluid_str = f"{v_fluid:<15.4f}" if v_fluid > 0 else f"{'[No Flow/Zero]':<15}"
        print(f"{cycle_time:<15} | {global_min:<11} | {dia:<15.3f} | {den:<10.3f} | {ws:<18.4f} | {v_fluid_str}")
        
    last_cycle_time = cycle_length - 1
    global_last_min = cycle_end_min
    if last_cycle_time % 10 != 0:
        dia = history["diameter"][global_last_min]
        den = history["density"][global_last_min]
        ws = history["wavespeed"][global_last_min]
        v_fluid = history["avg_fluid_vel"][global_last_min]
        
        if v_fluid > max_vel_in_cycle:
            max_vel_in_cycle = v_fluid
            
        v_fluid_str = f"{v_fluid:<15.4f}" if v_fluid > 0 else f"{'[No Flow/Zero]':<15}"
        print(f"{last_cycle_time:<15} | {global_last_min:<11} | {dia:<15.3f} | {den:<10.3f} | {ws:<18.4f} | {v_fluid_str} (End)")
    
    # Output the requested maximum tracking summary block at the bottom
    print("-" * 101)
    max_summary_str = f"{max_vel_in_cycle:.4f} µm/min" if max_vel_in_cycle > 0 else "[No Flow/Zero]"
    print(f"{'MAX VELOCITY IN CYCLE:':<78} {max_summary_str}")
    print("=" * 101)
