#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 26 10:35:47 2026
@author: maggie
"""
# DOT PRODUCT ANALYSIS: RFP (Mechanical) vs FRET/CFP (Biochemical)

from glob import glob
import os
import csv
import matplotlib.pyplot as plt
import numpy as np
from skimage.io import imread
from openpiv import pyprocess, validation, filters
from skimage.filters import gaussian, threshold_local, threshold_otsu
from skimage.draw import polygon
from skimage.util import img_as_float

# Declaring because OpenPiv says so :p
np.int = int

#paths
search_path_RFP = "/Volumes/qiongy-data/Users/Lucy/Data/04292026FreshEnergyInterphaseCycling/Pos10/img_*_7-RFP-T_000.tif"
search_path_FRET = "/Volumes/qiongy-data/Users/Maggie/OLD 04-29-26 PIV results/ImageJdata/POS10 ratio data/POS 10 ratio color/IMD_2.2_00*.tif"

rfp_paths = sorted(glob(search_path_RFP))
fret_paths = sorted(glob(search_path_FRET))

#checking that they are lined up 
num_frames = min(len(rfp_paths), len(fret_paths))
if num_frames < 2:
    raise ValueError(f"Insufficient matching frames found! RFP: {len(rfp_paths)}, FRET: {len(fret_paths)}")

rfp_paths = rfp_paths[:num_frames]
fret_paths = fret_paths[:num_frames]
total_pairs = num_frames - 1

print(f"Synchronized pipelines: Processing {total_pairs} paired frames.")

#setting up the PIV variables
dt = 3.0  # 3 minutes between frames
winsize = 32 
searchsize = 32  
overlap = 8 

vmin = 0.0
vmax = 5.0

# Create temporary holding space for saved visualization frames (wicked meme mention)
temp_dir = "combined_temp_frames"
os.makedirs(temp_dir, exist_ok=True)

#relativent tracking lists
time_list = []
piv_data_rows = []
frame_files = []

#the is for the final image (all from OPENPIV tut)
fig, ax = plt.subplots(figsize=(8, 8))
dummy_Q = ax.quiver([0], [0], [0], [0], [0], cmap="plasma", scale=115, clim=(vmin, vmax))
dummy_Q2 = ax.quiver([0], [0], [0], [0], [0], cmap="viridis", scale=115, clim=(vmin, vmax))
cb = fig.colorbar(dummy_Q, ax=ax, orientation='horizontal', pad=0.08)
cb2 = fig.colorbar(dummy_Q2, ax=ax, orientation='horizontal', pad=0.08)
cb.set_label('RFP Velocity Magnitude')
cb2.set_label('FRET/CFP Velocity Magnitude')
fig.tight_layout()

#okay this loop is gonna go through both folders RFP and FRET
print(f"Processing {total_pairs} frame pairs...")

#this is the same loop at algo_mask code
for i in range(total_pairs):
    print(f"Processing frame pair {i + 1} of {total_pairs}...")
    raw_fret_curr = imread(fret_paths[i])
    raw_fret_next = imread(fret_paths[i + 1])
    if raw_fret_curr.ndim == 3:
        fret_curr = img_as_float(raw_fret_curr[:, :, 0])
        fret_next = img_as_float(raw_fret_next[:, :, 0])
    else:
        fret_curr = img_as_float(raw_fret_curr)
        fret_next = img_as_float(raw_fret_next)
    fret_curr = np.clip(np.nan_to_num(fret_curr), 0.0, 2.0)
    fret_next = np.clip(np.nan_to_num(fret_next), 0.0, 2.0)
    #reading and adding next frame
    rfp_curr = img_as_float(imread(rfp_paths[i]))
    rfp_next = img_as_float(imread(rfp_paths[i + 1]))
    #denoising
    rfp_curr_denoised = gaussian(rfp_curr, sigma=1.0)
    rfp_next_denoised = gaussian(rfp_next, sigma=1.0)
    fret_curr_denoised = gaussian(fret_curr, sigma=1.0)
    fret_next_denoised = gaussian(fret_next, sigma=1.0)
    # OTSU stuff (google for more info)
    otsu_rfp_c = threshold_otsu(rfp_curr_denoised)
    otsu_rfp_n = threshold_otsu(rfp_next_denoised)
    otsu_fret_c = threshold_otsu(fret_curr_denoised)
    otsu_fret_n = threshold_otsu(fret_next_denoised)
    #this ^ and this (below) are pulled from the algo_mask code just now 2x for both sets of data
    rfp_c_img = (np.where(rfp_curr_denoised > otsu_rfp_c, rfp_curr, 0.0) * 255).astype(np.uint8)
    rfp_n_img = (np.where(rfp_next_denoised > otsu_rfp_n, rfp_next, 0.0) * 255).astype(np.uint8)
    fret_c_img = (np.where(fret_curr_denoised > otsu_fret_c, fret_curr, 0.0) * 255).astype(np.uint8)
    fret_n_img = (np.where(fret_next_denoised > otsu_fret_n, fret_next, 0.0) * 255).astype(np.uint8)
    #PIV for RFP
    u_rfp, v_rfp, s2n_rfp = pyprocess.extended_search_area_piv(
        rfp_c_img.astype(np.int16), rfp_n_img.astype(np.int16),
        window_size=winsize, overlap=overlap, dt=dt, search_area_size=searchsize, sig2noise_method='peak2peak')
    #PIV for FRET
    u_fret, v_fret, s2n_fret = pyprocess.extended_search_area_piv(
        fret_c_img.astype(np.int16), fret_n_img.astype(np.int16),
        window_size=winsize, overlap=overlap, dt=dt, search_area_size=searchsize, sig2noise_method='peak2peak')
    #this is to make sure the frames are lined up, coord geo too
    x, y = pyprocess.get_coordinates(image_size=rfp_c_img.shape, search_area_size=searchsize, overlap=overlap)
    y_indices, x_indices = y.astype(int), x.astype(int)
    #get rid of outliers
    flags_rfp = validation.sig2noise_val(u_rfp, v_rfp, s2n_rfp, threshold=1.0015)
    flags_fret = validation.sig2noise_val(u_fret, v_fret, s2n_fret, threshold=1.0015)
    u_rfp, v_rfp = filters.replace_outliers(u_rfp, v_rfp, flags_rfp, method='localmean', max_iter=10, kernel_size=2)
    u_fret, v_fret = filters.replace_outliers(u_fret, v_fret, flags_fret, method='localmean', max_iter=10, kernel_size=2)
    u_rfp, v_rfp = u_rfp.astype(float), v_rfp.astype(float)
    u_fret, v_fret = u_fret.astype(float), v_fret.astype(float)
    mask_layer = (rfp_curr_denoised <= otsu_rfp_c)
    img_poly1, img_poly2, img_poly3, img_poly4 = np.zeros_like(rfp_c_img, dtype=bool), np.zeros_like(rfp_c_img, dtype=bool), np.zeros_like(rfp_c_img, dtype=bool), np.zeros_like(rfp_c_img, dtype=bool)
    #hard setting a specific spot to look at 
    rr, cc = polygon([0, 0, 1150, 1150], 
                     [0, 400, 400, 0])
    qq, jj = polygon([0, 0, 560, 560], 
                     [0, 1150, 1150, 0])
    pp, kk = polygon([0, 0, 1150, 1150], 
                     [730, 1150, 1150, 730])
    mm, nn = polygon([900, 900, 1150, 1150],
                     [0, 1150, 1150, 0])
    img_poly1[rr, cc] = True
    img_poly2[qq, jj] = True
    img_poly3[pp, kk] = True
    img_poly4[mm,nn] = True
    #combine the algo and the hard mask
    grid_mask = (mask_layer[y_indices, x_indices] | img_poly1[y_indices, x_indices] | img_poly2[y_indices, x_indices] | img_poly3[y_indices, x_indices] | img_poly4[y_indices, x_indices])
    #defining both gigamasks 
    GIGAMASK_RFP = grid_mask | np.isnan(u_rfp) | np.isnan(v_rfp)
    GIGAMASK_FRET = grid_mask | np.isnan(u_fret) | np.isnan(v_fret)
    masked_u_rfp = np.ma.masked_array(u_rfp, mask=GIGAMASK_RFP)
    masked_v_rfp = np.ma.masked_array(v_rfp, mask=GIGAMASK_RFP)
    masked_u_fret = np.ma.masked_array(u_fret, mask=GIGAMASK_FRET)
    masked_v_fret = np.ma.masked_array(v_fret, mask=GIGAMASK_FRET)
    magnitude_rfp = np.ma.sqrt(masked_u_rfp**2 + masked_v_rfp**2)
    magnitude_fret = np.ma.sqrt(masked_u_fret**2 + masked_v_fret**2)
    
    
    
    #Okay this is where the dot product comes in
    #we are checking the angle between the RFP and FRET/CFP channels vectors in the same spot
    #mulitpling the arrays and then adding them, the dotproduct 
    dotprod = (masked_u_rfp * masked_u_fret) + (masked_v_rfp * masked_v_fret)
    #defining each of their magnitudes 
    norm_RFP = np.ma.sqrt(masked_u_rfp**2 + masked_v_rfp**2)
    norm_FRET = np.ma.sqrt(masked_u_fret**2 + masked_v_fret**2)
    #so theta = cos^-1(dot / magntitude of the 2 vectors)
    denom = norm_RFP * norm_FRET
    #adding limit
    cos_theta = np.clip(dotprod / denom, -1.0, 1.0)
    #converting it to degrees and taking the cos^-1
    theta_field_deg = np.degrees(np.ma.arccos(cos_theta))
    valid_thetas = theta_field_deg.compressed()
    current_time = dt * i
    for theta_val in valid_thetas:
        piv_data_rows.append([current_time, theta_val])
        
    #end of dot product stuff 
        
    time_list.append(dt * i)
    ax.clear()
    ax.imshow(rfp_c_img, alpha=0.5, cmap='gray')
    Q = ax.quiver(x, y, masked_u_rfp, -masked_v_rfp, magnitude_rfp, cmap="plasma", scale=115, width=0.005)
    Q = ax.quiver(x, y, masked_u_fret, -masked_v_fret, magnitude_fret, cmap="viridis", scale=115, width=0.005)
    Q.set_clim(vmin, vmax)
    
    ax.invert_yaxis() 
    ax.set_xlim(0, rfp_c_img.shape[1])
    ax.set_ylim(rfp_c_img.shape[0], 0)
    ax.set_title(f"Velocity Field Cross-Analysis | Minute: {int(dt * i)}")
  
    frame_path = os.path.join(temp_dir, f"frame_{i:03d}.png")
    plt.savefig(frame_path, dpi=100)
    frame_files.append(frame_path)

print("\nAll frames evaluated and rendered successfully.")
plt.close()
csv_filename = "/Users/maggie/Desktop/Time and theta_POS10.csv"
print(f"Saving metric alignments down to {csv_filename}...")
with open(csv_filename, mode='w', newline='') as f:
    writer = csv.writer(f)
    writer.writerow(["Time(min)", "Theta(Degrees)"])
    writer.writerows(piv_data_rows)
print("CSV analytics saved successfully!")
print("Cleaning up temporary local assets...")
for f in frame_files:
    os.remove(f)
os.rmdir(temp_dir)
print("Process finished safely. End of Code.")