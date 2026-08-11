#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 11 14:38:04 2026

@author: maggie
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

#YOUR DATA HERE!
x_data = np.array([
    572, 574, 576, 578, 580, 582, 584, 586, 588, 590,
    592, 594, 596, 598, 600, 602, 604, 606, 608, 610,
    612, 614, 616, 618, 620, 622, 624, 626, 628, 630,
    632, 634, 636, 638, 640, 642, 644, 646, 648, 650
])

y_data = np.array([
      0,   0,   0,   0,   0,   5,  20,  40,  75, 120,
    130, 140, 145, 150, 160, 150, 165, 165, 160, 165,
    165, 165, 165, 170, 170, 160, 165, 160, 165, 165,
    160, 165, 165, 170, 170, 180, 180, 180, 185, 175
])

# sigmodal function
def sigmoid(x, bottom, top, x0, k):
    return bottom + (top - bottom) / (1 + np.exp(-k * (x - x0)))

#initialzing
p0 = [min(y_data), max(y_data), np.median(x_data), 0.1]

#creating our fit.
popt, _ = curve_fit(sigmoid, x_data, y_data, p0=p0)
bottom, top, x0, k = popt

#create our threshold.
threshold = 0.99
y_plateau = bottom + threshold * (top - bottom)

#solving for where our 99% threshold is met. 
x_plateau = x0 - (1 / k) * np.log((top - bottom) / (y_plateau - bottom) - 1)

#printing out results. 
print("--- Fit Results ---")
print(f"Lower Asymptote (Bottom Plateau Y): {bottom:.2f}")
print(f"Upper Asymptote (Max Plateau Y):    {top:.2f}")
print(f"Inflection Point (Midpoint X):      {x0:.2f}")
print(f"\n--- 99% Plateau Point ---")
print(f"X (where plateau starts): {x_plateau:.2f}")
print(f"Y (at plateau threshold): {y_plateau:.2f}")

#presentation stuff 
x_fit = np.linspace(min(x_data), max(x_data), 300)
y_fit = sigmoid(x_fit, *popt)
plt.figure(figsize=(9, 5))
plt.scatter(x_data, y_data, color="black", label="Data Points", zorder=3)
plt.plot(x_fit, y_fit, color="red", linewidth=2, label="Sigmoidal Fit")
plt.axhline(top, color="blue", linestyle="--", alpha=0.7, label=f"Max Asymptote Y ≈ {top:.1f}")
plt.plot(x_plateau, y_plateau, "go", markersize=8, label=f"99% Plateau (X={x_plateau:.1f}, Y={y_plateau:.1f})")
plt.xlabel("X")
plt.ylabel("Y")
plt.title("Sigmoidal Fit and Plateau Detection")
plt.legend()
plt.grid(True, linestyle=":", alpha=0.6)
plt.tight_layout()
plt.show()