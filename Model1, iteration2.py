#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 15 16:19:43 2026

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
T = [74, 56, 60, 66, 60, 58, 64, 64, 118, 110, 108]
t_T = [122, 178, 238, 304, 364, 422, 486, 550, 668, 778, 886]

k = random.randint(5, 15)

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
    "wavespeed": []
}

effective_diameter = 0.001 #microns, simple tubulin

current_cycle = 0
cycle_time = 0


#iterate through all 1200 minutes
for i in range(0, t, 1):
    
    #determining which cycle we are in dependning on how much time has elapsed
    if current_cycle < len(t_T) and i == t_T[current_cycle]:
        current_cycle += 1
        cycle_time = 0
    
    #cycles 0-2
    if 0 <= current_cycle < 3:
        if cycle_time == 0:
            effective_diameter = 0.001
            polyrate = 0
        elif 1 <= cycle_time <= 35:
            polyrate = (10.583 * cycle_time - 68.421) * 2.6
            if polyrate > 0:
                effective_diameter = polyrate * random.uniform(.9, 1.1)
        elif 35 < cycle_time <= 90:
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
        elif 24 < cycle_time <= 90:
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
        elif 38 < cycle_time <= 90:
            polyrate = (0.0296 * cycle_time + 212.39) * 2.6
            if polyrate > 0:
                effective_diameter = polyrate * random.uniform(.9, 1.1)
                
    #calculations
    current_volume = volume(effective_diameter)
    l_density = local_density(current_volume)
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

cycle_starts = [0] + t_T[:-1]
cycle_ends = [val - 1 for val in t_T]

total_cycles = len(T)

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
    print(f"{'Time in Cycle':<15} | {'Global Min':<11} | {'Diameter (µm)':<15} | {'Density':<10} | {'Wavespeed (µm/min)':<18}")
    print("-" * 79)
    for cycle_time in range(0, cycle_length, 10):
        global_min = cycle_start_min + cycle_time
        dia = history["diameter"][global_min]
        den = history["density"][global_min]
        ws = history["wavespeed"][global_min]
        print(f"{cycle_time:<15} | {global_min:<11} | {dia:<15.3f} | {den:<10.3f} | {ws:<18.4f}")
    last_cycle_time = cycle_length - 1
    global_last_min = cycle_end_min
    if last_cycle_time % 10 != 0:
        dia = history["diameter"][global_last_min]
        den = history["density"][global_last_min]
        ws = history["wavespeed"][global_last_min]
        print(f"{last_cycle_time:<15} | {global_last_min:<11} | {dia:<15.3f} | {den:<10.3f} | {ws:<18.4f} (End)")
    print("-" * 79)
    