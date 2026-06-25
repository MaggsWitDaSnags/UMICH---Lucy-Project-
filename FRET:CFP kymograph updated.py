#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 22 15:51:14 2026

@author: maggie
"""
#HARD MASK VERSION
#comment out the stuff you need

#library calls
import csv
import os
from glob import glob
import matplotlib.pyplot as plt
import numpy as np
from skimage.io import imread
from skimage.util import img_as_float
from skimage.filters import gaussian, threshold_otsu
from skimage.draw import polygon

#function declaration for processing the data and masking
#note that these y/x mins/maxs can be altered depending what you would like your window to look like
def process_and_mask_frame(raw_img, ymin=300, ymax=1100, xmin=0, xmax=400):
    
    #we must check if the incoming images are in B&W, otherwise we must convert them
    if raw_img.ndim == 3:
        img_float = img_as_float(raw_img[:, :, 0])
    else:
        img_float = img_as_float(raw_img)
    
    #decarling the img as a float in order to alter it
    img_float = np.nan_to_num(img_float, nan=0.0, posinf=0.0, neginf=0.0)
    img_float = np.clip(img_float, 0.0, 2.0)

    #a lot of this was taken directly from a tutorial
    #overall this chunk preforms a gaussian calculation to find the spot with the largest change and denoises the rest 
    img_denoised = gaussian(img_float, sigma=1.0)
    #this sets the threshold for our masking, anything below that point, we leave out i.e. darkness around the sample
    otsu_val = threshold_otsu(img_denoised)
    #and setting that threshold true for background noise
    intensity_mask = (img_denoised <= otsu_val) 
    
    #YOUR sandbox to play around with creating 'hard masks' and not auto ones
    #This should be used for specific isolation
    #google how to properly use the polygon function for more details
    #the more masks you need, the more mask decalrations you will need
    poly_mask1 = np.zeros_like(img_float, dtype=bool)
    poly_mask2 = np.zeros_like(img_float, dtype=bool)
    #poly_mask3 = np.zeros_like(img_float, dtype=bool)
    
    #parameters for the mask, these are like coordinates
    rr1, cc1 = polygon([0, 346, 346, 0],[0, 0, 1150, 1150])
    rr2, cc2 = polygon([0, 1150, 1150, 0],[330, 330, 1150, 1150])
    #rr3, cc3 = polygon([0, 1150, 1150, 0], [0, 0, 292, 292])
    
    #set those areas to true in order for them to block out analysis 
    poly_mask1[rr1, cc1] = True
    poly_mask2[rr2, cc2] = True
    #poly_mask3[rr3, cc3] = True
    
    #the whole mask = gigamask
    #composed of the algo and the hard masks 
    #using OR operator and opposed to &
    gigamask = intensity_mask | poly_mask1 | poly_mask2 #| poly_mask3
    
    #declaring our mask area as where the gigamask defines on our image
    masked_img = np.where(gigamask, 0.0, img_float)
    
    #returning that image cropped at our min/max values
    return masked_img[ymin:ymax, xmin:xmax]

#path delcarations 
#mac is funny with path names, google if unsure
search_path = "/Users/maggie/Desktop/POS12ratiodata/color/IMD_1.8_00*.tif"
#taking in all of our images at the path and subracting off frame that are not neccesary 
image_paths = sorted(glob(search_path))[:-185]

#making sure we got the images 
#note if you are using this code and are not MAGGIE, you can always throw in print check throughout the code to catch were an error might be occuring 
if len(image_paths) < 2:
    raise ValueError("Double check your path configuration!")
print(f"Successfully tracked {len(image_paths)} ratio frames for mitotic wave analysis.")

#setting our time step 
dt = 3.0  # 3 minutes between frames
#since this code looks at changes in pairs, we must declare our total number of pairs to parse through
total_pairs = len(image_paths) - 1

#reading first image and checking if it is B&W, if not changing 
first_raw = imread(image_paths[0])
first_img = img_as_float(first_raw[:, :, 0]) if first_raw.ndim == 3 else img_as_float(first_raw)
#decalring our shape of the image, neccesary for kymograph data 
height, width = first_img.shape

#declaring our array where we will store realvent data 
all_diffs_2d = [] 
time_list = [] 
global_activity = [] 

#ANOTHA PRINT CHECK
print(f"Processing {total_pairs} frame pairs with active GIGAMASK filters...")

#now, our image under analysis has the proper parameters: color change and masked
curr_frame_border = process_and_mask_frame(imread(image_paths[0]))

#loop where caluclations come from
for i in range(total_pairs):
    
    #counting the time
    current_time = i * dt
    time_list.append(current_time)
    
    # Read ahead to next frame and process/mask it
    raw_next = imread(image_paths[i + 1])
    next_frame_border = process_and_mask_frame(raw_next)
    
    # Absolute difference calculation (now noise-free!)
    #this is where the global change data comes from 
    diff_frame = np.abs(next_frame_border - curr_frame_border)
    #add this to our array 
    all_diffs_2d.append(diff_frame)
    #taking that mean difference and adding it to our global activity array 
    global_activity.append(np.mean(diff_frame))
    
    #Setting up loop for the next frame v frame analysis 
    curr_frame_border = next_frame_border

#creating an array in the current array 
all_diffs_2d = np.array(all_diffs_2d)

#Generate Kymographs
# taking the median values (yeilded the clearest data, also at suggestion of lab)
#taking the axis = 1, x-axis and transposing the data 
kymo_x_wave = np.median(all_diffs_2d, axis=1).T 
#similar done here for y-axis 
kymo_y_wave = np.median(all_diffs_2d, axis=2).T  


#directional slope calc 
with np.errstate(divide='ignore', invalid='ignore'):
    min_len = min(kymo_y_wave.shape[0], kymo_x_wave.shape[0])
    delta_ratio_kymo = kymo_y_wave[:min_len, :] / (kymo_x_wave[:min_len, :] + 1e-6)

#Locate Wave Source Origin
#this is still a work in progress 
#maping the overall variance in the all_diffs_2d array
spatial_variance_map = np.var(all_diffs_2d, axis=0)
#declaring that teh osurce must come from the point of max change overtime
wave_y_source, wave_x_source = np.unravel_index(np.argmax(spatial_variance_map), spatial_variance_map.shape)

#We must adjust these variables if we are cropping the image, so do so here
ymin, xmin = 300, 0
global_wave_x = wave_x_source + xmin
global_wave_y = wave_y_source + ymin

#print check!
print(f"\nCalculated Mitotic Wave Origin: X={global_wave_x}, Y={global_wave_y}")

#EXPORTING DATA TO CSV --> need to workshop...
csv_filename = "wave_source_analytics.csv"
print(f"Saving compiled spatial metrics to {csv_filename}...")
with open(csv_filename, mode='w', newline='') as f:
    writer = csv.writer(f)
    writer.writerow(["Time (min)", "Global Wave Activity Metric", "Mean Kymo X Profile", "Mean Kymo Y Profile"])
    writer.writerows(zip(time_list, global_activity, np.mean(kymo_x_wave, axis=0), np.mean(kymo_y_wave, axis=0)))
    

#PLOT 1.X-Kymograph
plt.figure(figsize=(10, 4))
plt.imshow(kymo_x_wave, cmap='inferno', aspect='auto', extent=[0, (total_pairs * dt), 0, kymo_x_wave.shape[0]])
plt.colorbar(label='Wave Activity ($\Delta$ Ratio)')
plt.xlabel('Time (minutes)')
plt.ylabel('Cropped X Position (pixels)')
plt.title('Mitotic Wave X-Kymograph (Median Projection with GIGAMASK)')
plt.tight_layout()
plt.savefig('wave_kymograph_xpos12FRET_CFP.png', dpi=150)
plt.show()

#PLOT 2.Directional Slope Kymograph
plt.figure(figsize=(10, 4))
plt.imshow(delta_ratio_kymo, cmap='twilight', aspect='auto', vmin=0, vmax=3, extent=[0, (total_pairs * dt), 0, min_len])
plt.colorbar(label='Propagation Ratio (Y/X)')
plt.xlabel('Time (minutes)')
plt.ylabel('Spatial Index')
plt.title('Wave Directional Dynamics: $\Delta$Y / $\Delta$X')
plt.tight_layout()
plt.savefig('wave_ratio_slopepos12FRET_CFP.png', dpi=150)
plt.show()

#PLOT 3.Wave Origin Map, work in prog 
plt.figure(figsize=(6, 6))
plt.imshow(first_img, cmap='gray')
#plt.scatter(global_wave_x, global_wave_y, marker='*', color='red', s=150, label='Wave Origin')
plt.imshow(spatial_variance_map, cmap='jet', alpha=0.4,) #extent=[xmin, xmin + spatial_variance_map.shape[1], ymin + spatial_variance_map.shape[0], ymin])
plt.title('2D Masked Wave Activity Map')
plt.colorbar(label='Temporal Variance')
plt.legend()
plt.tight_layout()
plt.savefig('wave_source_origin_mappos12FRET_CFP.png', dpi=150)
plt.show()

#final print check, runs when code has completed. 
print("\nAnalysis complete! Masked metrics and plots generated successfully.")