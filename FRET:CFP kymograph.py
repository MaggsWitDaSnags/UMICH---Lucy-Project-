#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jun 21 17:43:15 2026

@author: maggie
"""
# FRET/CFP CODE!!!!!!!!
# USE THIS ONE IF YOU DON't NEED TO HARD MASK
import csv
import os
from glob import glob
import matplotlib.pyplot as plt
import numpy as np
from skimage.io import imread
from skimage.util import img_as_float

#path declaring, google if issues arise, mac is annoying with paths
search_path = "/Users/maggie/Desktop/POS14ratiodata/color/IMD_2.1_00*.tif"
#cutting off how many frmaes we take in if they are not useful
image_paths = sorted(glob(search_path))[:-185]

#double checks that path is correct, a print check
if len(image_paths) < 2:
    raise ValueError("Double check your path configuration!")
print(f"Successfully tracked {len(image_paths)} ratio frames for mitotic wave analysis.")

dt = 3.0  # 3 minutes between frames
#since we are looking at images in pairs and the change between them, we need to declare the number of pairs
total_pairs = len(image_paths) - 1

# Read first frame to get dimensions
#we must convert the image to B&W if the data is in color, 3D-->2D
first_raw = imread(image_paths[0])
if first_raw.ndim == 3:
    first_img = img_as_float(first_raw[:, :, 0])
else:
    first_img = img_as_float(first_raw)
#hard setting the shape
height, width = first_img.shape

#pre-declaring the kymographs, width = x-axis, height = y-axis
kymo_x_wave_T = np.zeros((total_pairs, width))  
kymo_y_wave_T = np.zeros((total_pairs, height)) 

#for calculating running variance without storing all frames: Var(X) = E[X^2] - (E[X])^2
sum_diff = np.zeros((height, width), dtype=np.float64)
sum_diff_sq = np.zeros((height, width), dtype=np.float64)

#declaring the neccesary arrays to store our data in 
time_list = []
global_activity = []

print(f"Processing {total_pairs} frame pairs to capture dynamic wavefront shifts...")

# Loop where values are calculated 
for i in range(total_pairs):
    #added this in bc it was running slowwww. only loading in 10 frames at a time, not 900!
    if (i + 1) % 10 == 0 or i == 0 or i == total_pairs - 1:
        print(f"Processing frame {i + 1} of {total_pairs}...")
    
    #declaring our time and adding it to our time array 
    current_time = i * dt
    time_list.append(current_time)
    
    #decalring both our current and latter image
    raw_curr = imread(image_paths[i])
    raw_next = imread(image_paths[i + 1])
    
    # Changing into 2d array, incoming images are color
    if raw_curr.ndim == 3:
        curr_frame = img_as_float(raw_curr[:, :, 0])
        next_frame = img_as_float(raw_next[:, :, 0])
    else:
        curr_frame = img_as_float(raw_curr)
        next_frame = img_as_float(raw_next)
      
    #finding and declaring our values that may be 0 or infity to be finite, division errors prevenot   
    curr_frame = np.nan_to_num(curr_frame, nan=0.0, posinf=0.0, neginf=0.0)
    next_frame = np.nan_to_num(next_frame, nan=0.0, posinf=0.0, neginf=0.0)
    #hard setting our range of values acceptable 0 --> 2
    curr_frame = np.clip(curr_frame, 0.0, 2.0)
    next_frame = np.clip(next_frame, 0.0, 2.0) 
    
    #caluclating the aboslute change between frames
    diff_frame = np.abs(next_frame - curr_frame)
    
    #Fill kymograph components for this specific time frame
    #taking median change, this approach was suggested by the lab 
    kymo_x_wave_T[i, :] = np.median(diff_frame, axis=0) #declaring our axis
    kymo_y_wave_T[i, :] = np.median(diff_frame, axis=1) #declaring our axis
    
    #the mean net change between frames is our global activity
    global_activity.append(np.mean(diff_frame))
    
    #total sum our our differences = roughly the average sqaured 
    #we need this for caluclations out side of the loop
    #calcualtions that look at big picture 
    sum_diff += diff_frame
    sum_diff_sq += diff_frame ** 2

#finalize Kymographs (transpose back to your original expected shapes: space vs time)
kymo_x_wave = kymo_x_wave_T.T
kymo_y_wave = kymo_y_wave_T.T

#caluclating the mean differnce over the total number of pairs
mean_diff = sum_diff / total_pairs
#the squared differences
mean_diff_sq = sum_diff_sq / total_pairs
#the spatial difference calulation the mean - the mean diff^2
spatial_variance_map = mean_diff_sq - (mean_diff ** 2)
#clean up any potential tiny negative floats due to floating-point rounding errors
spatial_variance_map = np.clip(spatial_variance_map, 0, None)

#directional slope calc 
with np.errstate(divide='ignore', invalid='ignore'):
    min_len = min(kymo_y_wave.shape[0], kymo_x_wave.shape[0])
    delta_ratio_kymo = kymo_y_wave[:min_len, :] / (kymo_x_wave[:min_len, :] + 1e-6)

#trying to find most staurated point, theoretically the area of orgin 
#WORK IN PROGRESS
wave_y_source, wave_x_source = np.unravel_index(np.argmax(spatial_variance_map), spatial_variance_map.shape)
print(f"\nCalculated Mitotic Wave Origin: X={wave_x_source}, Y={wave_y_source}")

#creating csv
csv_filename = "wave_source_analytics.csv"
print(f"Saving compiled spatial metrics to {csv_filename}...")
with open(csv_filename, mode='w', newline='') as f:
    writer = csv.writer(f)
    writer.writerow(["Time (min)", "Global Wave Activity Metric", "kymoxwave", "kymoywave"])
    writer.writerows(zip(time_list, global_activity, kymo_x_wave, kymo_y_wave))
    
#x-based kymograph 
plt.figure(figsize=(10, 4))
plt.imshow(kymo_x_wave, cmap='inferno', aspect='auto', extent=[0, (total_pairs * dt), 0, width])
plt.colorbar(label='Wave Activity ($\Delta$ Ratio)')
plt.xlabel('Time (minutes)')
plt.ylabel('X Position (pixels)')
plt.title('Mitotic Wave X-Kymograph (Median Projection)')
plt.tight_layout()
plt.savefig('wave_kymograph_xpos13FRET_CFP.png', dpi=150)
plt.show()

#slope-based kymograph
plt.figure(figsize=(10, 4))
plt.imshow(delta_ratio_kymo, cmap='twilight', aspect='auto', vmin=0, vmax=3, extent=[0, (total_pairs * dt), 0, min_len])
plt.colorbar(label='Propagation Ratio (Y/X)')
plt.xlabel('Time (minutes)')
plt.ylabel('Spatial Index')
plt.title('Wave Directional Dynamics: $\Delta$Y / $\Delta$X')
plt.tight_layout()
plt.savefig('wave_ratio_slopepos13FRET_CFP.png', dpi=150)
plt.show()

#wave source assumption map
#WORK IN PROGRESS...dont trust yet
plt.figure(figsize=(6, 6))
plt.imshow(first_img, cmap='gray')
#plt.scatter(wave_x_source, wave_y_source, )#marker='*', color='red')
plt.imshow(spatial_variance_map, cmap='jet', alpha=0.4)
plt.title('2D Wave Activity Map')
plt.colorbar(label='Temporal Variance')
plt.tight_layout()
plt.savefig('wave_source_origin_mappos13FRET_CFP.png', dpi=150)
plt.show()

#print showing the end of teh sim 
print("\nAnalysis complete! All files and plots are generated and saved.")