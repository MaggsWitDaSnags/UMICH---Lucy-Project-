#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 28 14:09:19 2026

@author: maggie
"""
#Model 2 --> Iteration 2

"""
Hypothesis:
    
    The mechanics of tubulin polyermization and depolyermization intitated by the CDK1 wave
    differ in compartment forming samples and non-compartment forming samples. 
    
    In compartment forming samples the tubulin form compartments or clusters that are dense
    in nature and can trap activation or inhibitory protiens/signals. The force that is 
    extered onto the cytoplasm once these cluster are depolyermized by the cdk1 signal, 
    which is sent out as a wave, impose a force larger enough to create a homogeneous mixture
    within the cytoplasm:
        
        Heterogeneous (compartments, tubulin segmentation) --> Homogeneous (well mixed)
        
    The mixing of the indhibitory and activating proteins in the compartment forming samples 
    maintain a relatively constant period until arrest (this is due to resource depletion).
    As opposed to non-compartment forming samples, the periodicty slowly increase with time, 
    most likley do to both resource depletion and the lack of uniform mixing the the sample.
    
    Therefore, it can be hypothesized that the compartments/clusters act as spherical springs
    in the sample that affect the cytoplasmic flow and the effective diffusion of the active
    signaling proteins/molecules. This downstream then also affects transportation.
    
    Overall, the periocidty, wavespeed and phaselag between mechanical actvity and mitotic 
    activity can be characterized in this model as a result of the mechanical forces 
    impossed by the tubulin during depolyermization. 
    
    The phase shift and periodicty data will be used a validation and compared to experimental 
    results. 
    
    
Parameters: 
    Polyrate --> measured from microscope data
    Diamater --> the result of polyrate 
    Radi --> calculated from diameter 
    Delta Radii --> caluclated via radii 
    Growth Rate --> calculated by the change in radii over a given time
    Stored Mechanical Energy --> Measured by the change in radii compared to the rest radii 
                                 *This can be changed into depend on growth rate later 
    Mechanical Force --> Measured by the change in energy over the change in radii 
    'Mechanical Release' --> defines the new period starting and the old period ending 
    Mixed Fraction --> a function of the released mechanical force post depolyermzation of tubulin
    Mixing Threshold --> This threshold a function of the initial mixed fraction post 
                         depolyermzation. I think later period can be just defined by the mixed 
                         fraction, this step may be unneccesary 
    Mixedness --> this is a function to help determine diffusion, I think this should be altered...?
    Diffusion --> This is a function of the mixedness which should decrease as the tubulin polyermize
                  due to the clusters creating barriers? not sure about this one just yet... needs to
                  be workshopped. 
    Wave Speed --> not implemented yet but, v = 2sqrtD??? OR a function of transport
    Transport --> This is a function of diffusion and of cytoplasmic flow, this eq should stay 
    Delta Phase Shift --> This is a function of transport, from V1. I think this  simple relation
                          should stay since it was 'validated' in the previous model, V1.
    
Assumptions:
    - Does directionallity of the force matter?
    - I think later we should assume that the compartments decrease in size overtime
      no longer using the polyrate numbers?
    - A lot of the basline variables...
    - seperating the flow and diffusion to net into transport
      The flow is mechanically induced and the diffusion is chemically induced?
    - Treating the cluster as a spring
    
"""

#global libraries
import numpy as np
import matplotlib.pyplot as plt


#global variables 
t = 1200 #minutes, total time under microscope 
dt = 1 #timestep 
k_eff = 1e-6  #Mechanical Stiffness (N/microns)
rest_radii = 10.0 #microns, can change later 
max_force_for_pure_homogeneous_mixture = 0.2 #Newton, can adjust later 
baseline_drag = 1.0
baseline_flow = 1.0
baseline_diffusion = 1.0
beta = 1.0
gamma = 1.0
theta = 1.0
mixed_threshold = 1.0
constant = 1.0


#Functions
def growth_rate(delta_rad):
    
    if delta_rad == 0:
        growth_rate = 0
    else:
        growth_rate = max(delta_rad,0) / dt
    
    return growth_rate

def calculated_mechanical_energy_stored(radii_parameter): #should I put growth rate here?
    
    stored_mech_energy = (0.5 * (k_eff) * ((radii_parameter - rest_radii)**2)) 
    
    return stored_mech_energy

def stored_energy_calc(arr):
    
    return max(arr) - min(arr)


def calculated_local_flow(force_parameter):
    
    local_cyto_flow = max(baseline_flow, baseline_flow + ((beta * force_parameter) - (gamma * baseline_drag)))
    
    return local_cyto_flow

def updated_diffuision_calc(mixing):
    
    diffusion = baseline_diffusion*(1+theta*mixing) # I think in this case it should decrease?
    # like it will decrease as the compartments form so maybe a whole different calcualtoin for mixing?
    return diffusion

def calculated_transport(diffusion_parameter, flowparameter):
    
    transport = (diffusion_parameter / baseline_diffusion) * (flowparameter / baseline_flow)
    
    return transport

def calculated_delta_phase_shift(transport_parameter):
    
    if transport_parameter == 0:
        phaseshift = np.pi
    else:
        phaseshift = np.pi / transport_parameter
    
    return phaseshift


#arrays to store my data 
minute = []
curr_diameter = []
curr_radius = []
delta_radii = []
growth_r8 = []
stored_mechanical_E = []
delta_stored_energy = []
mechanical_force = []
T = [60,]
cyto_local_flow = []
mixing_val = []
stored_mechanical_E_cycle = []
diffusion = []
transport = []
delta_phase_shift = []

#initalizing
current_cycle = 0
cycle_time = 0
Period = 60 #mins
sample_mixed_fraction = 0.70 #starts off pretty homogeneous 

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
            max_rad = 266.5 + 20 #microns
            if cycle_time == 0:
                effective_diameter = 0.008
                effective_radius = effective_diameter/2
                polyrate = 0
            elif 1 <= cycle_time <= 42:
                polyrate = (7.6703 * cycle_time - 22.956) * 2.6
                effective_diameter = polyrate
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
     
    delta_radii.append(delta_radii_carrier)
        
    #growth rate stuff
    growth_rate_carrier = growth_rate(delta_radii_carrier)
    growth_r8.append(growth_rate_carrier)
    
    #stored mechanical enegry stuff
    stored_mech_energy_carrier = calculated_mechanical_energy_stored(effective_radius)
    stored_mechanical_E.append(stored_mech_energy_carrier)
    stored_mechanical_E_cycle.append(stored_mech_energy_carrier)
    
    if i==0:
        delta_energy_carrier = 0
    else:
        delta_energy_carrier = stored_mechanical_E[i] - stored_mechanical_E[i-1]
    
    delta_stored_energy.append(delta_energy_carrier)
    
    released_energy = stored_energy_calc(stored_mechanical_E_cycle)
    
    mechanical_force_carrier = released_energy * constant
    mechanical_force.append(mechanical_force_carrier)
    sample_mixed_fraction = abs(mechanical_force_carrier / max_force_for_pure_homogeneous_mixture)
   
    #Mechanical release --> the follow occur after the mechanical release
    if cycle_time >= Period:

        mechanical_force_carrier = released_energy * constant 
        mechanical_force[i] = mechanical_force_carrier
        sample_mixed_fraction = abs(mechanical_force_carrier / max_force_for_pure_homogeneous_mixture)
    
        current_cycle += 1
        cycle_time = 0
        Period = 60 + 15*(1-sample_mixed_fraction)
        T.append(Period)
        stored_mechanical_E_cycle.clear()
        
    mixing_val.append(sample_mixed_fraction)
    
    
    #local flow stuff 
    local_cyto_flow_carrier = calculated_local_flow(released_energy)
    cyto_local_flow.append(local_cyto_flow_carrier)
    
    #diffusion stuff 
    diffusion_carrier = updated_diffuision_calc(sample_mixed_fraction)#need to work this calc
    diffusion.append(diffusion_carrier)
    
    #transport stuff
    transport_carrier = calculated_transport(diffusion_carrier, local_cyto_flow_carrier)
    transport.append(transport_carrier)
    
    #phase stuff 
    delta_phase_shift_carrier = calculated_delta_phase_shift(transport_carrier)
    delta_phase_shift.append(delta_phase_shift_carrier)
    
    

#printing stuff 
for j in range(280,360,dt):
    print("MINUTE: ", minute[j], 
          "Radii: ", curr_radius[j],
          "Diameter: ", curr_diameter[j],
          "Change in radii: ", delta_radii[j],
          "Stored Mechanical Energy: ", stored_mechanical_E[j],
          "Cluster Growth Velocity: ", growth_r8[j],
          "Local Cytoplasmic Flow: ", cyto_local_flow[j],
          "Diffusion: " ,diffusion[j],
          "Transport:", transport[j],
          "Delta Phase Shift", delta_phase_shift[j])
    
for m in range(0,len(T),dt): 
     print(T[m])
     

    
plt.plot(minute, curr_diameter , label = 'diameter')
plt.legend() 
for x in range(60,1201, 60):
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

plt.plot(minute, stored_mechanical_E , color = 'red', label = 'stored mechanical Energy')
plt.legend()  
for x in range(60,1201,60):
    plt.axvline(x, linestyle='--', alpha=0.3, color='black') 
plt.show()

plt.plot(minute, mechanical_force , color = 'blue', label = 'Mechanical force') #gotta fix this so that it only peaks at the cdk1 wave
plt.legend()  
for x in range(60,1201,60):
    plt.axvline(x, linestyle='--', alpha=0.3, color='black') 
plt.show()

#plt.plot(minute, T, color = 'blue', label = 'Period')                                
#plt.legend()  
#for x in range(60,1201,60):
    #plt.axvline(x, linestyle='--', alpha=0.3, color='black')                             
#plt.show()

plt.plot(minute, cyto_local_flow, color = 'darkorange', label = 'Local Cytoplasmic Flow' )
plt.legend() 
for x in range(60,1201,60):
    plt.axvline(x, linestyle='--', alpha=0.3, color='black')  
plt.show()

plt.plot(minute, mixing_val, color = 'red', label = 'Mixing Threshold')
plt.legend()  
for x in range(60,1201,60):
    plt.axvline(x, linestyle='--', alpha=0.3, color='black') 
plt.show()

plt.plot(minute, diffusion, color = 'cyan', label = 'Diffusion')
plt.legend()  
for x in range(60,1201,60):
    plt.axvline(x, linestyle='--', alpha=0.3, color='black') 
plt.show()
    
plt.plot(minute, transport, color = 'purple', label = 'transport coeffeient')
plt.legend()  
for x in range(60,1201,60):
    plt.axvline(x, linestyle='--', alpha=0.3, color='black') 
plt.show()

plt.plot(minute, delta_phase_shift, color = 'gray', label = 'change in phase shift')
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
