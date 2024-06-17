from helper_calculations.landscape import Configuration
from helper_calculations.fisher_information import build_map_FIMS, plot_uncertainty_ellipse
from helper_calculations.make_map import make_basic_seismic_map
from helper_calculations.make_target_map import add_targets_localized
import numpy as np
import pprint as pp
from pathlib import Path
sensor_comm_ratio = 1.5
import matplotlib.pyplot as plt
import matplotlib

# -------- DESCRIPTION -----------------
# This script was used to plot fig 2 in the 2024 Sensors
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
def add_sensor(ax, position, radius, sensor_type, solt_type, init):
    if sensor_type == "radial":
        if init == "initial":
            color_pick = "darkblue"
            fill = True
        elif init == "final":
            color_pick = "darkblue"
            fill = False
    elif sensor_type == "bearing":
        color_pick = "#ff7f00"
        fill = True
    
    if solt_type == "analytical":
        line_style = "-"
        hatch = ""
    elif solt_type == "numerical":
        line_style = "--"
        hatch = '///////'
    elif solt_type == "numerical_LOS":
        line_style = "--"
        hatch = '+++++'

    # Add the sensing radius
    if sensor_type == "bearing":
        circle = plt.Circle((position[0], position[1]), radius, color = color_pick, linewidth = 2, fill=False, linestyle = line_style, lw = 1.5)
        ax.add_patch(circle)
        circle2 = plt.Circle((position[0], position[1]), 1.5*radius, color = color_pick, linewidth = 2, fill=False, linestyle = "--", lw = 1)
        ax.add_patch(circle2)

    # Add the sensor
    #if sensor_type == "radial":
    circle = plt.Circle((position[0], position[1]), 4, color = color_pick, linewidth = 2, fill=fill, linestyle = '-', lw = 1, hatch = hatch)
    ax.add_patch(circle) 

    return ax

## CREATE FIG 1/2 - DIRECT VS ANALYTICAL 
sensor_locs_initial = [52, 35, 52, 51]
sensor_rad = [15, 15, 15]
num_sensors = len(sensor_rad)
meas_type = ["bearing", "radial", "radial"]
targets = [49, 42]
localized = [1] # corresponds to number of target
sensor_locs_DIRECT = [40.129629629629626, 45.833333333333336] #w/o LOS placement
sensor_locs_analytical = [41.01432225726775, 45.2342691911711]

#LOS VALUES DIRECT 
LOS_flag = 1 #1 = on, 0 off

# CREATE THE FIMS FOR ELLIPSES
# Analytical
sensor_locs_analytical = [52, 35, 52, 51, sensor_locs_analytical[0], sensor_locs_analytical[1]]
inputs = localized, targets, sensor_locs_analytical, sensor_rad, meas_type, terrain, LOS_flag
(FIMs_analytical, _) = build_map_FIMS(inputs)

# DIRECT - No LOS
sensor_locs_FIM_DIRECT = [52, 35, 52, 51, sensor_locs_DIRECT[0], sensor_locs_DIRECT[1]]
inputs = localized, targets, sensor_locs_FIM_DIRECT, sensor_rad, meas_type, terrain, LOS_flag
(FIMs_DIRECT, _) = build_map_FIMS(inputs)

# DIRECT - w/ LOS
sensor_locs_DIRECT_LOS = [55.611111111111114, 40.129629629629626] #w/ LOS placement
sensor_locs_FIM_DIRECT = [52, 35, 52, 51, sensor_locs_DIRECT_LOS[0], sensor_locs_DIRECT_LOS[1]]
inputs = localized, targets, sensor_locs_FIM_DIRECT, sensor_rad, meas_type, terrain, LOS_flag
(FIMs_DIRECT_LOS, _) = build_map_FIMS(inputs)

# ----------------------------------

# Create a plot with targets and sensors
sensor_locs_analytical = [6*i for i in sensor_locs_analytical]
targets = [6*i for i in targets]
sensor_rad = [6*i for i in sensor_rad]
sensor_locs_DIRECT = [6*i for i in sensor_locs_DIRECT]
sensor_locs_DIRECT_LOS = [6*i for i in sensor_locs_DIRECT_LOS]

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

# Add the analytical sensors to the map
sensor_type = ["bearing","radial","radial"]
inits = ["initial", "initial", "final"]
for k in range(len(sensor_locs_analytical)//2):
    ax = add_sensor(ax, [sensor_locs_analytical[0+2*k], sensor_locs_analytical[1+2*k]], sensor_rad[k], sensor_type[k], "analytical", inits[k])

# Add the DIRECT sensor to the map
for k in range(len(sensor_locs_DIRECT)//2):
    ax = add_sensor(ax, [sensor_locs_DIRECT[0+2*k], sensor_locs_DIRECT[1+2*k]], sensor_rad[k], "radial", "numerical", "final")

# Add the DIRECT-LOS sensor to the map
for k in range(len(sensor_locs_DIRECT_LOS)//2):
    ax = add_sensor(ax, [sensor_locs_DIRECT_LOS[0+2*k], sensor_locs_DIRECT_LOS[1+2*k]], sensor_rad[k], "radial", "numerical_LOS", "final")


# Add the target
ax = add_targets_localized(ax, targets, localized)

# Add the uncertainty ellipses
# Analytical
ax = plot_uncertainty_ellipse(ax, FIMs_analytical[0], targets, 1.52, 6, "black", "analytical")
print(FIMs_analytical)

# DIRECT - Naive
ax = plot_uncertainty_ellipse(ax, FIMs_DIRECT[0], targets, 1.52, 6, "magenta", "numerical")
print(FIMs_DIRECT)

# DIRECT - LOS
ax = plot_uncertainty_ellipse(ax, FIMs_DIRECT_LOS[0], targets, 1.52, 6, "black", "numerical")
print(FIMs_DIRECT_LOS)


# Formatting
#ax.set_ylabel('Longitude'), ax.set_xlabel('Latitude')
ax.set_title('')
#ax.legend(["_bearing rad", "_Bearing Sensor", "_Distance Sensor", "_", "Exact Placement", "DIRECT Placement", 
#           "_target", "Exact $1\sigma$", "DIRECT $1\sigma$"])
ax.legend(["_1", "_1", "_1", "_1", "Naive-Exact Placement", "Naive-Numerical Placement", 
          "Env-Aware Placement", "_target", "Naive-Exact 1$\sigma$", "Naive-Numerical 1$\sigma$", "Env-Aware 1$\sigma$", "_1"])
plt.xlim([150, 400]), plt.ylim([150, 350])
plt.gca().invert_yaxis()

plt.show()