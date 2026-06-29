#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 25 10:33:59 2026

@author: maggie
"""
#This is the belnded code of both the PIV analysis and FRET/CFP (kymograph & global activity) data
#library call
import os
import csv
from glob import glob
import numpy as np
import matplotlib.pyplot as plt
import imageio.v3 as iio
from skimage.io import imread
from skimage.util import img_as_float
from skimage.filters import gaussian, threshold_otsu
from openpiv import pyprocess, validation, filters

#delcaring paths for both of the folders containing relavent data 
# Kept local on your Mac Desktop
rfp_path_pattern = "/Volumes/qiongy-data/Users/Lucy/Data/04292026FreshEnergyInterphaseCycling/Pos14/img_*_7-RFP-T_000.tif"

# UPDATED: Connected directly to your 100 GB network server folder
fret_path_pattern = "/Volumes/qiongy-data/Users/Maggie/OLD 04-29-26 PIV results/ImageJdata/POS14ratiodata/color/IMD_2.1_00*.tif"

#can add [:-###] if you want to remove x number of frames from teh back end
rfp_paths = sorted(glob(rfp_path_pattern))
fret_paths = sorted(glob(fret_path_pattern))

#lining up all of the images
#print check to ensure that the frames are lined up and there is enough data  
num_frames = min(len(rfp_paths), len(fret_paths))
if num_frames < 2:
    raise ValueError(f"Insufficient matching frames found! RFP: {len(rfp_paths)}, FRET: {len(fret_paths)}")

#creating variables for loops
rfp_paths = rfp_paths[:num_frames]
fret_paths = fret_paths[:num_frames]
total_pairs = num_frames - 1

#print check
print(f"Synchronized pipelines: Processing {total_pairs} paired frames ({num_frames} total frames).")

dt = 3.0  # 3 minutes between frames

#size of the window we are looking at 
winsize = 32 
searchsize = 32  
overlap = 8 

#hard setting the velocity scale, mostly for asetetics 
vmin, vmax = 0.0, 5.0  # Vector velocity scale

#Loading in the fret images to determine the image height and width
# FIXED: Added indexing [0] so imread looks at the first file string instead of the list object
first_fret_raw = imread(fret_paths[0])
fret_height, fret_width = first_fret_raw.shape[:2]

#this is for efficentcy, I suck at spelling holy 
temp_dir = "combined_temp_frames"
os.makedirs(temp_dir, exist_ok=True)

#declaring our lists where the data will be stored
time_list = []
avg_speed_list = []
max_speed_list = []
global_fret_activity = []

#setting up data containers for the kymographs 
kymo_x_wave_T = np.zeros((total_pairs, fret_width))  
kymo_y_wave_T = np.zeros((total_pairs, fret_height)) 

#presentation stuff
fig, ax = plt.subplots(figsize=(8, 8))
#this us for the arrows, see other code for reference
dummy_Q = ax.quiver([0], [0], [0], [0], [0], cmap="plasma", scale=115, clim=(vmin, vmax))
cb = fig.colorbar(dummy_Q, ax=ax, orientation='horizontal', pad=0.08)
cb.set_label('Mechanical Velocity Magnitude')

#another data container
frame_files = []

#actual frame processing code
#currently, on my mac, this is running for about 10~15 minutes 
#the logic in this loop is almost identical to the other code 
for i in range(total_pairs):
    if (i + 1) % 10 == 0 or i == 0 or i == total_pairs - 1:
        print(f"Processing paired frame {i + 1} of {total_pairs}...")

    current_time = i * dt
    time_list.append(current_time)

    curr_rfp = imread(rfp_paths[i])
    next_rfp = imread(rfp_paths[i + 1])
    
    raw_fret_curr = imread(fret_paths[i])
    raw_fret_next = imread(fret_paths[i + 1])

    #converting color to  B&W
    if raw_fret_curr.ndim == 3:
        fret_curr = img_as_float(raw_fret_curr[:, :, 0])
        fret_next = img_as_float(raw_fret_next[:, :, 0])
    else:
        fret_curr = img_as_float(raw_fret_curr)
        fret_next = img_as_float(raw_fret_next)

    fret_curr = np.clip(np.nan_to_num(fret_curr), 0.0, 2.0)
    fret_next = np.clip(np.nan_to_num(fret_next), 0.0, 2.0)

    #fret data 
    fret_diff = np.abs(fret_next - fret_curr)
    global_fret_activity.append(np.mean(fret_diff))
    kymo_x_wave_T[i, :] = np.median(fret_diff, axis=0)
    kymo_y_wave_T[i, :] = np.median(fret_diff, axis=1)

    #RFP data 
    rfp_curr_float = img_as_float(curr_rfp)
    rfp_next_float = img_as_float(next_rfp)

    rfp_curr_blur = gaussian(rfp_curr_float, sigma=1.0)
    rfp_next_blur = gaussian(rfp_next_float, sigma=1.0)

    otsu_curr = threshold_otsu(rfp_curr_blur)
    otsu_next = threshold_otsu(rfp_next_blur)

    rfp_curr_prep = (np.where(rfp_curr_blur > otsu_curr, rfp_curr_float, 0.0) * 255).astype(np.int16)
    rfp_next_prep = (np.where(rfp_next_blur > otsu_next, rfp_next_float, 0.0) * 255).astype(np.int16)

    mask_layer = (rfp_curr_blur <= otsu_curr)

    u, v, sig2noise = pyprocess.extended_search_area_piv(
        rfp_curr_prep,
        rfp_next_prep,
        window_size=winsize,
        overlap=overlap,
        dt=dt,
        search_area_size=searchsize,
        sig2noise_method='peak2peak',
    )

    x, y = pyprocess.get_coordinates(image_size=rfp_curr_prep.shape, search_area_size=searchsize, overlap=overlap)
    flags = validation.sig2noise_val(u, v, sig2noise, threshold=1.0015)
    u, v = filters.replace_outliers(u, v, flags, method='localmean', max_iter=10, kernel_size=2)

    u, v = u.astype(float), v.astype(float)
    grid_mask = mask_layer[y.astype(int), x.astype(int)]
    
    GIGAMASK = grid_mask | np.isnan(u) | np.isnan(v)
    masked_u = np.ma.masked_array(u, mask=GIGAMASK)
    masked_v = np.ma.masked_array(v, mask=GIGAMASK)
    magnitude = np.ma.sqrt(masked_u**2 + masked_v**2)

    avg_speed_list.append(np.nanmean(magnitude))
    max_speed_list.append(np.nanmax(magnitude))

   # HEAT MAP CODE 
   #ax.clear()
    #okay this plots the FRET/CFP data as a heat map 
    #heatmap = ax.imshow(
        #fret_diff, 
        #cmap='inferno', 
        #alpha=0.9, 
        #vmin=0.0, 
        #vmax=0.5  # <-- play around with this depending on incoming data 
    #)
    
    #overlaying the RFP arrows 
    #Q = ax.quiver(
        #x, y, 
        #masked_u, -masked_v,
        #magnitude, 
        #cmap="plasma", 
        #scale=115, # <-- play around with this also depending on imcoming data 
        #width=0.005
    #)
    #Q.set_clim(vmin, vmax)
    
    #more presentation stuff, making it all stay aligned 
    #keeps the coords aligned 
    #ax.set_xlim(0, fret_width)
    #ax.set_ylim(fret_height, 0) 
    #ax.set_title(f"PIV Vectors over FRET Heatmap | Minute: {int(current_time)}")

    #frame_path = os.path.join(temp_dir, f"overlay_{i:03d}.png")
    #plt.savefig(frame_path, dpi=100)
    #frame_files.append(frame_path) #adding the frames to the master data container
    
    #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    #if you just wanna do the raw FRET/CFP images:
    #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ax.clear()
    ax.imshow(raw_fret_curr)
    Q = ax.quiver(
        x, y, 
        masked_u, -masked_v,
        magnitude,
        cmap="plasma", 
        scale=115, 
        width=0.005
    )
    Q.set_clim(vmin, vmax)
    ax.set_xlim(0, fret_width)
    ax.set_ylim(fret_height, 0) 
    ax.set_title(f"PIV Vectors over Raw FRET/CFP Image | Minute: {int(current_time)}")
    frame_path = os.path.join(temp_dir, f"overlay_{i:03d}.png")
    plt.savefig(frame_path, dpi=100)
    frame_files.append(frame_path)
plt.close()

#data output stuff 
print("\nStitching synchronized overlay movie...")
images = [iio.imread(f) for f in frame_files]
gif_output_path = "/Users/maggie/Desktop/mechanical_biochemical_overlayRAWMAPPOS14.gif"
iio.imwrite(gif_output_path, images, plugin="pillow", duration=250, loop=0)

#csv output, can always add more data to this in the future depending on what we want to see 
csv_filename = "/Users/maggie/Desktop/synchronized_wave_analyticsPOS14.csv"
print(f"Saving compiled spatial metrics to {csv_filename}...")
with open(csv_filename, mode='w', newline='') as f:
    writer = csv.writer(f)
    writer.writerow(["Time(min)", "AvgMechanicalSpeed", "MaxMechanicalSpeed", "GlobalBiochemicalFRETDelta"])
    writer.writerows(zip(time_list, avg_speed_list, max_speed_list, global_fret_activity))

#more cleanin 
print("Cleaning up temp files...")
for f in frame_files:
    os.remove(f)
os.rmdir(temp_dir)

#indcated end of code
#another print check 
print("\nPipeline execution complete! Overlay saved as 'mechanical_biochemical_overlay.gif'")
