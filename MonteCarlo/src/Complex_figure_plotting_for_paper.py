from helper_calculations.landscape import Configuration, Node
from helper_calculations.sensors import Sensor
import numpy as np
import pprint as pp
import os
import shutil
from pathlib import Path
from scipy.optimize import minimize
from scipy import optimize
import time
from helper_calculations.make_map import make_basic_seismic_map
from helper_calculations.make_target_map import add_targets, add_targets_localized
sensor_comm_ratio = 1.5
import matplotlib.pyplot as plt
from helper_calculations.fisher_information import plot_uncertainty_ellipse, build_FIM, build_map_FIMS
import matplotlib
import random

# -------- DESCRIPTION -----------------
# This script was used to plot fig 3 in the 2024 Sensors
# conference paper we submitted. 
# --------------------------------------

# -------- Make LaTeX kinda fonts -----------------
matplotlib.rcParams['mathtext.fontset'] = 'custom'
matplotlib.rcParams['mathtext.rm'] = 'Bitstream Vera Sans'
matplotlib.rcParams['mathtext.it'] = 'Bitstream Vera Sans:italic'
matplotlib.rcParams['mathtext.bf'] = 'Bitstream Vera Sans:bold'

matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'


# IMPORT TERRAIN
# Import the 100x100 where we do our calculations
terrain_height = 100 #796
terrain_width = 100 #1002
terrain = Configuration(terrain_width, terrain_height)
my_path = Path(__file__).parent / "../data/terrain/GIS_terrain_resize.csv"
terrain.load_from_csv(my_path)


# Def a function to add the circles and sqrs to the plot
def add_sensor(ax, position, radius, sensor_type, solt_type, sens_rad_plot):
    
    # Colors and plotting for different cnds
    if solt_type == "Naive":
        hatch = ''
        fill = False
        
    elif solt_type == "initial":
        hatch = '///////'
        fill = True
        
    elif solt_type == "Env":
        hatch = '///////'
        fill = False
        color_pick = "#377eb8" #blue

    if sensor_type == "bearing":
        color_pick = "#ff7f00" #orange
    elif sensor_type == "radial" or "radius":
        color_pick = "darkblue" #blue

    # Place the sensors on the plot
    # Make the initial sensors squares and the final sensors circles
    # Keep colors corresponding to the sensor type
    plot_rad = 3.5
    #if solt_type == "initial":
    #    w, h = 2*plot_rad, 2*plot_rad
    #    rect = plt.Rectangle((position[0]-w/2, position[1]-h/2), w, h, color = color_pick, linewidth = 2, fill=fill, linestyle = '-', lw = 1.5, hatch = hatch)
    #    ax. add_patch(rect)

    if solt_type == "Naive" or "Env" or "initial":
        # add a black circle for legend purposes
        circle = plt.Circle((position[0], position[1]), plot_rad, color = "black", fill=fill, linestyle = '-', lw = 1, hatch = hatch)
        ax.add_patch(circle)

        circle = plt.Circle((position[0], position[1]), plot_rad, color = color_pick, linewidth = 2, fill=fill, linestyle = '-', lw = 1.5, hatch = hatch)
        ax.add_patch(circle)

    # Add sensing rad if necessary
    if sens_rad_plot == 1:
        # Add the sens radius
        circle = plt.Circle((position[0], position[1]), radius, color = color_pick, linewidth = 2, fill=False, linestyle = '-', lw = 1)
        ax.add_patch(circle)
        # Add the comm radius
        circle = plt.Circle((position[0], position[1]), radius*1.5, color = color_pick, linewidth = 2, fill=False, linestyle = '--', lw = 1)
        ax.add_patch(circle)

    return ax

## INITIAL SENSOR LIST
sensor_rad = [25, 25, 25]
sensor_type = ["seismic","acoustic","seismic"]
num_sensors = len(sensor_type)
sensor_comm_ratio = 2 # ratio of sensor communication to sensing radius 
meas_type = ["radius", "bearing", "radius"] #radius or bearing
sensor_locs = []
sens1 = [29, 44] 
sens2 = [33, 29] 
sens3 = [59, 35] 
#sens4 = [44, 29] 
sensor_locs.extend(sens1)
sensor_locs.extend(sens2)
sensor_locs.extend(sens3)
#sensor_locs.extend(sens4)

## TARGET LIST
targets = []
tar1 = [37, 39]
tar2 = [42, 32] 
tar3 = [47, 39] 
tar4 = [41, 45] 
targets.extend(tar1)
targets.extend(tar2)
targets.extend(tar3)
targets.extend(tar4)

localized = [1,1,1,1]
LOS_flag = 1 #1 = on, 0 off


## NAIVE FINAL PLACEMENT
# For 3 added sensor solts - 6 units away

sensor_rad_new = [15, 15, 20]
sensor_type_new = ["seismic", "acoustic", "seismic"]
num_sensors_new = len(sensor_type_new)
sensor_comm_ratio_new = [2] # ratio of sensor communication to sensing radius 
meas_type_new = ["radial", "bearing","radial"]
sensor_rad.extend(sensor_rad_new)
sensor_type.extend(sensor_type_new)
meas_type.extend(meas_type_new)

sensor_locs_naive = [49.9, 45.833333333333336, 48.27777777777778, 32.38888888888889, 46.64814814814815, 51.944444444444436]
sensor_locs_env = [30.759259259259263, 38.5, 50, 32.8, 48.27777777777778, 23.833333333333336]


# For 4 added sensor solts
'''
sensor_rad_new = [15, 20, 15, 20]
sensor_type_new = ["seismic", "seismic", "acoustic", "acoustic"]
num_sensors_new = len(sensor_type_new)
sensor_comm_ratio_new = [2] # ratio of sensor communication to sensing radius 
meas_type_new = ["radial", "radial","bearing","bearing"]
sensor_rad.extend(sensor_rad_new)
sensor_type.extend(sensor_type_new)
meas_type.extend(meas_type_new)

sensor_locs_naive = [32.38, 51.12, 48.27, 51.13, 48.27, 25.87, 47.87, 51.12]
sensor_locs_env = [27.5, 40.94, 34.83, 21.38, 44.61, 58.05, 48.27, 23.83]


# For 3 added sensor solts - 8 units away
sensor_rad_new = [15, 15, 20]
sensor_type_new = ["seismic", "acoustic", "seismic"]
num_sensors_new = len(sensor_type_new)
sensor_comm_ratio_new = [2] # ratio of sensor communication to sensing radius 
meas_type_new = ["radial", "bearing","radial"]
sensor_rad.extend(sensor_rad_new)
sensor_type.extend(sensor_type_new)
meas_type.extend(meas_type_new)

sensor_locs_naive = [32.38, 51.13, 48.27, 25.87, 48.27, 51.13]
sensor_locs_env = [27.5, 40.94, 44.61, 58.05, 48.27, 23.83]
'''

# CREATE THE FIMS FOR ELLIPSES
# Naive
sensor_locs_naive_input = sensor_locs.copy()
sensor_locs_naive_input.extend(sensor_locs_naive)
print(sensor_locs_naive_input)
inputs = localized, targets, sensor_locs_naive_input, sensor_rad, meas_type, terrain, LOS_flag
(FIMs_naive, det_naive) = build_map_FIMS(inputs)
print("Naive det:", det_naive)
print(FIMs_naive)

# Env
sensor_locs_env_input = sensor_locs.copy()
sensor_locs_env_input.extend(sensor_locs_env)
print(sensor_locs_env_input)
inputs = localized, targets, sensor_locs_env_input, sensor_rad, meas_type, terrain, LOS_flag
(FIMs_env, det_env) = build_map_FIMS(inputs)
print("Env det:", det_env)
print(FIMs_env)

# ----------------------------------

# Create a plot with targets and sensors
sensor_locs_naive = [6*i for i in sensor_locs_naive]
sensor_locs_env = [6*i for i in sensor_locs_env]
sensor_locs = [6*i for i in sensor_locs]
targets = [6*i for i in targets]
sensor_rad = [6*i for i in sensor_rad]


# IMPORT TERRAIN
# Import the 600 x 600 to actually make the plot
terrain_height = 600 #796
terrain_width = 600 #1002
terrain = Configuration(terrain_width, terrain_height)
my_path = Path(__file__).parent / "../data/terrain/GIS_terrain_cropped.csv"
terrain.load_from_csv(my_path)

# Creates the basic colormap
fig = plt.figure()
_map_init = make_basic_seismic_map(0, [], [], [], sensor_comm_ratio, [])
ax = terrain.plot_grid(_map_init)

# Add the initial sensors to the map
sens_plot_rad_list = [0, 0, 0] # use 1 to plot the sensor radius and the comm radius
for k in range(len(sensor_locs)//2):
    # ax, position, radius, sensor_type, solt_type, init
    ax = add_sensor(ax, [sensor_locs[0+2*k], sensor_locs[1+2*k]], sensor_rad[k], meas_type[k], "initial", sens_plot_rad_list[k])

# Add the newly placed sensors
sens_plot_rad_list = [0, 0, 0, 0] # use 1 to plot the sensor radius and the comm radius
for k in range(len(sensor_locs_naive)//2):
    # ax, position, radius, sensor_type, solt_type, init
    ax = add_sensor(ax, [sensor_locs_naive[0+2*k], sensor_locs_naive[1+2*k]], sensor_rad[k], meas_type_new[k], "Naive", sens_plot_rad_list[k])
    ax = add_sensor(ax, [sensor_locs_env[0+2*k], sensor_locs_env[1+2*k]], sensor_rad[k], meas_type_new[k], "Env", sens_plot_rad_list[k])

# Add annotations
# Do the fixed sensors first
#plt.text(sensor_locs[0]+0, sensor_locs[1]-6, "$R_{1,F}$", fontsize = 'large', fontweight='bold')  
#plt.text(sensor_locs[2]+4, sensor_locs[3]-6, "$B_{F}$", fontsize = 'large', fontweight='bold')  
#plt.text(sensor_locs[4]+2, sensor_locs[5]-6, "$R_{2,F}$", fontsize = 'large', fontweight='bold') 

# 3 sensor placement map

# Do the Naive sensors 
#plt.text(sensor_locs_naive[0]+3, sensor_locs_naive[1]-5, "$R_{1,N}$", fontsize = 'x-large', fontweight='bold')  
#plt.text(sensor_locs_naive[2]+2, sensor_locs_naive[3]-4, "$B_{1,N}$", fontsize = 'x-large', fontweight='bold')  
#plt.text(sensor_locs_naive[4]+3, sensor_locs_naive[5]-3, "$R_{2,N}$", fontsize = 'x-large', fontweight='bold')  

# Do the Env sensors 
#plt.text(sensor_locs_env[0]+2, sensor_locs_env[1]+9, "$R_{1,E}$", fontsize = 'x-large', fontweight='bold')  
#plt.text(sensor_locs_env[2]+3, sensor_locs_env[3]+2, "$B_{1,E}$", fontsize = 'x-large', fontweight='bold')  
#plt.text(sensor_locs_env[4]+3, sensor_locs_env[5]-3, "$R_{2,E}$", fontsize = 'x-large', fontweight='bold')  


# 4 sensor placement map
'''
# Do the Naive sensors 
plt.text(sensor_locs_naive[0]+3, sensor_locs_naive[1]-5, "$R_{1,N}$", fontsize = 'x-large', fontweight='bold')    
plt.text(sensor_locs_naive[2]+3, sensor_locs_naive[3]-3, "$R_{2,N}$", fontsize = 'x-large', fontweight='bold')  
plt.text(sensor_locs_naive[4]+2, sensor_locs_naive[5]-4, "$B_{1,N}$", fontsize = 'x-large', fontweight='bold')
plt.text(sensor_locs_naive[6]+2, sensor_locs_naive[7]-4, "$B_{2,N}$", fontsize = 'x-large', fontweight='bold')

# Do the Env sensors 
plt.text(sensor_locs_env[0]+2, sensor_locs_env[1]+9, "$R_{1,E}$", fontsize = 'x-large', fontweight='bold')   
plt.text(sensor_locs_env[2]+3, sensor_locs_env[3]-3, "$R_{2,E}$", fontsize = 'x-large', fontweight='bold') 
plt.text(sensor_locs_env[4]+3, sensor_locs_env[5]+2, "$B_{1,E}$", fontsize = 'x-large', fontweight='bold') 
plt.text(sensor_locs_env[6]+3, sensor_locs_env[7]+2, "$B_{2,E}$", fontsize = 'x-large', fontweight='bold') 
'''

# Add the target
ax = add_targets_localized(ax, targets, localized)

# Add the uncertainty ellipses
for i, _ in enumerate(FIMs_env):
    target = [targets[0+2*i], targets[1+2*i]]
    # Naive
    ax = plot_uncertainty_ellipse(ax, FIMs_naive[i], target, 2.4, 6, "black", "numerical")
    # Env
    ax = plot_uncertainty_ellipse(ax, FIMs_env[i], target, 2.4, 6, "black", "analytical")
    

# Formatting
#ax.set_ylabel('Longitude'), ax.set_xlabel('Latitude')
ax.set_title('')
ax.legend(["_1", "Fixed Range Sensor", "_1", "Fixed Bearing Sensor", "_1", "_1", "Naive Placed Sensor", "_1", "Env-Aware Placed Sensor", 
           "_1", "_1", "_1", "_1", "_1", "_1", "_1", "_1", "_1", "_1", "_1", "_1", 
           "_1", "Naive 95% CI", "_1", "_1", "_1", "_1", "Env-Aware 95% CI"], loc='lower right')

#ax.legend(["Fixed Range", "Fixed Bearing", "_3", "Naive Placement", "Env.-Aware Placement", 
#           "_6", "_7", "_8", "_9", "_10", "_11", "_12", "_13", "Naive 95% CI", "Env.-Aware 95% CI", "_16", "_17"], loc='lower right')
plt.xlim([125, 425]), plt.ylim([100, 350])
plt.gca().invert_yaxis()

# Annotate 3 configs
#ax.annotate('Initial Bearing', xy=(198, 174), xytext=(185, 119), arrowprops={'arrowstyle':'-|>'})
#ax.annotate('Naively-Placed Bearing', xy=(292, 192), xytext=(308, 158), arrowprops={'arrowstyle':'-|>'})
#ax.annotate('Env.-Placed Radial', xy=(290, 142), xytext=(320, 123), arrowprops={'arrowstyle':'-|>'})

plt.show()