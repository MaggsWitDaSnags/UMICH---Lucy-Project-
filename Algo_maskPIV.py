#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 15 16:37:44 2026

@author: maggie
"""
#Use for RFP only!!!!!! you can use for FRET/CFP but results are iffy...
#nor do they really make any sense bc they are a biochemical phenomena, not a mechanical one

#SOURCED BY: https://openpiv.readthedocs.io/en/latest/src/masking.html
#global library 
from glob import glob
import os
import imageio.v3 as iio
import matplotlib.pyplot as plt
import numpy as np
from skimage.io import imread
from openpiv import pyprocess, validation, filters, preprocess
from skimage.filters import gaussian, threshold_local
from skimage.util import img_as_float, img_as_ubyte
from skimage.filters import gaussian, threshold_otsu
from skimage.draw import polygon

#declaring cause the code got anrgy
np.int = int

#path set up 
#change 'POS9' to your folder of interest
#change 'RFP' to filter of intrest 
search_path = "/Users/maggie/Desktop/POS14/img_*_7-RFP-T_000.tif"
image_paths = sorted(glob(search_path))

#PRINT CHECK
#makes sure all files are being read
if len(image_paths) < 2:
    raise ValueError(
        f"Found {len(image_paths)} images at path: {search_path}\n"
        "Double check your path!"
    )

#ANOTHA PRINT CHECK
print(f"Successfully tracked {len(image_paths)} frames for processing.")

#reading first image
first_img = imread(image_paths[0])
#creating a temperary holding space for the files (wicked meme mention)
temp_dir = "temp_frames"
os.makedirs(temp_dir, exist_ok=True)

#hard setting our scale
vmin=0.0
vmax=5.0

#plotting the scale bar
fig, ax = plt.subplots(figsize=(8, 8))
dummy_Q = ax.quiver([0], [0], [0], [0], [0], cmap="plasma", scale=115, clim=(vmin,vmax))
cb = fig.colorbar(dummy_Q, ax=ax, orientation='horizontal', pad=0.08)
cb.set_label('Velocity Magnitude')
fig.tight_layout()

#bin for loop, collects the images into an array
frame_files = []
total_pairs = len(image_paths) - 1

# Storage lists for overall analytics
time_list = []
avg_speed_list = []
max_speed_list = []
avg_move_x_list = []
avg_move_y_list = []

#print check ot ensure the code is running up to here 
#parsing through all of the images
print(f"Processing {total_pairs} frame pairs...")

#LOOP WHERE BIG STUFF HAPPENS
for i in range(total_pairs):
    #counter for the images being processes
    print(f"Processing frame {i + 1} of {total_pairs}...")
    
    #THIS FOLLOWING IS CONFUSING...
    #SOURCED BY NUMPY and SKIMAGE
    
    #layering the frames ontop of one another to measure pixel change
    curr_frame = imread(image_paths[i])
    next_frame = imread(image_paths[i + 1])
    
    #scaling our array to more standardized ranges (0-->1)
    curr_float = img_as_float(curr_frame)
    next_float = img_as_float(next_frame)
    
    #blurring inbetween bc we just need to see movement in general, inc sigma for more blur
    curr_denoised = gaussian(curr_float, sigma=1.0)
    next_denoised = gaussian(next_float, sigma=1.0)
    
    #changes threshold value per round of images, no global value bc that messes it up!
    curr_thresh = threshold_local(curr_denoised, block_size=35, method='gaussian', offset=0.05)
    next_thresh = threshold_local(next_denoised, block_size=35, method='gaussian', offset=0.05)
    
    # Otsu threshold --> see wiki page &skimage page for more details
    otsu_curr = threshold_otsu(curr_denoised)
    otsu_next = threshold_otsu(next_denoised)
    
    #catches errors, sometimes the value will be larger than 1.0, in that case this fixes that 
    curr_frame = (np.where(curr_denoised > otsu_curr, curr_float, 0.0) * 255).astype(np.uint8)
    next_frame = (np.where(next_denoised > otsu_next, next_float, 0.0) * 255).astype(np.uint8)
    
    #print(f"Min: {next_float.min()}, Max: {next_float.max()}")
    #^comment this back in if code magically stops...
    
    #creating the mask to eliminate noise
    #checking where that mask should fall, anything that is not denoised
    mask_layer = (curr_denoised <= otsu_curr)

    #how many pixels we want to evaluate at one time
    winsize = 32 
    searchsize = 32  
    overlap = 8 
    dt = 3.0 #3min between images 
    
    #This stuff is a bit confusuing
    #please reference openPIV, masking tutuorial for more info
    #basically analyzing both images and converts then to 16x16 pictcure
    #Notices overlap or lack their of and creates a 2D array for both the x&y direction
    #ALL OF THIS STUFF HAS TO DO WITH VECTOR CALCULATIONS REFER TO OPENPIV FOR MORE DETAILS
    
    
    u, v, sig2noise = pyprocess.extended_search_area_piv(
        curr_frame.astype(np.int16),
        next_frame.astype(np.int16),
        window_size=winsize,
        overlap=overlap,
        dt=dt,
        search_area_size=searchsize,
        sig2noise_method='peak2peak',
    )

    #''
    x, y = pyprocess.get_coordinates(
        image_size=curr_frame.shape,
        search_area_size=searchsize,
        overlap=overlap,
    )
    
   # CONTROLING CRAZY VECTORS
    flags = validation.sig2noise_val(u, v, sig2noise, threshold=1.0015)
    
    #manual masking, comment out when neccesary 
    #use sandbox to test coords
    #img = np.zeros_like(curr_frame, dtype=bool)
    #img2 = np.zeros_like(curr_frame, dtype=bool)
    #img3 = np.zeros_like(curr_frame, dtype=bool)
    
    #rr, cc = polygon(
       #[0, 346, 346, 0],       
        #[0, 0, 1150, 1150],     
    #)
    
    #qq, jj = polygon (
        #[1100, 1150, 1150, 1010], 
        #[0, 0, 1150, 1150],       
    #)

    #yy, xx = polygon(
        #[0, 1150, 1150, 0],      
        #[330, 330, 1150, 1150],        
    #)
    
    #img[rr, cc] = True
    #img2[qq, jj] = True
   # img3[yy, xx] = True         
    
    y_indices = y.astype(int)
    x_indices = x.astype(int)
    
    #combining all of our masks (GIGAMASK tehe)
    grid_mask = mask_layer[y_indices, x_indices] #( img[y_indices, x_indices] | img3[y_indices, x_indices] | mask_layer[y_indices, x_indices]) #| img3[y_indices, x_indices]

    u, v = filters.replace_outliers(
        u, v, flags, method='localmean', max_iter=10, kernel_size=2,
    )
    
    # Force float arrays so we can safely inject NaNs
    u = u.astype(float)
    v = v.astype(float)
    
    
    # HARD MASK FORCE: If the grid mask is True (background), set vector to NaN
    u[grid_mask] = np.nan
    v[grid_mask] = np.nan
    
    GIGAMASK = grid_mask | np.isnan(u) | np.isnan(v)
    
    masked_u = np.ma.masked_array(u, mask=GIGAMASK)
    masked_v = np.ma.masked_array(v, mask=GIGAMASK)
    
    #calculating the magnitiude of the vectors
    magnitude = np.ma.sqrt(masked_u**2 + masked_v**2)

    # Appending calculation metrics for current frame pair
    current_time = dt * i
    time_list.append(current_time)
    avg_speed_list.append(np.nanmean(magnitude))
    max_speed_list.append(np.nanmax(magnitude))
    avg_move_x_list.append(np.nanmean(masked_u))
    avg_move_y_list.append(np.nanmean(masked_v))

   #PRESENTATION
    #start fresh after loop!
    ax.clear()
    #placing microscope image underneath 
    ax.imshow(curr_frame, alpha=0.5, cmap='gray')
    #draws in arrows
    Q = ax.quiver(
        x, 
        y, 
        masked_u, 
        -masked_v, 
        magnitude, 
        cmap="plasma",
        scale=115, 
        width=0.005,
    )
    Q.set_clim(vmin, vmax)
    #sacles colorbar incase or weirdness
    #cb.update_normal(Q)
    #more presentation stuff
    ax.invert_yaxis() 
    ax.set_xlim(0, first_img.shape[1])
    ax.set_ylim(first_img.shape[0], 0)
    ax.set_title(f"Velocity Field: {search_path } Minute: {3 * i}")
  
    frame_path = os.path.join(temp_dir, f"frame_{i:03d}.png")
    plt.savefig(frame_path, dpi=100)
    frame_files.append(frame_path)

#print check!
print("\nAll individual vector frames rendered successfully.")
plt.close()
#stiching all of the new images together to form a gif
print("Stitching frames together...")
images = [iio.imread(f) for f in frame_files]
#lower duration for faster GIF
iio.imwrite("velocity_field_movie_RFPPOS14algo.gif", images, plugin="pillow", duration=250, loop=0)
import csv
csv_filename = "velocity_analyticsRFPPOS14algo.csv"
#PRINT CHECK
print(f"Saving data to {csv_filename}...")
# putting the list we decalred before into rows
with open(csv_filename, mode='w', newline='') as f:
    writer = csv.writer(f)
    # headers
    writer.writerow(["Time (min)", "Avg Speed", "Max Speed", "Avg Move X", "Avg Move Y"])
    writer.writerows(zip(time_list, avg_speed_list, max_speed_list, avg_move_x_list, avg_move_y_list))
#PRINT CHECK
print("CSV saved successfully!")
#housekeeping...not overloading computa with unneccesary files
print("Cleaning up temp files...")
for f in frame_files:
    os.remove(f)
os.rmdir(temp_dir)
#print check!
print("\nMovie saved successfully as 'velocity_field_movie_RFPPOS14algo.gif'!")