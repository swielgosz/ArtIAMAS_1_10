from landscape import Configuration, Node
from sensors import Sensor
import numpy as np
import pprint as pp
import os
import shutil
from pathlib import Path
from scipy.optimize import minimize
import matplotlib.pyplot as plt
import seaborn as sns
from fisher_information import build_FIM

'''
This script is recreating some of the figures in "optimal analysis of sensor-
target localization geometries" (2010) for just the bearing only measurements

It fixes 2 sensors and then calculates the det(FIM) at every point in the grid
to create a heatmap of "good" and bad placements. 

NOTE: to avoid quick decay of the FIM, change the FIM.py script to that ri is
always 1! Or scale them down or something. Else, these plots aren't illuminating

'''

# TERRAIN INPUTS
terrain_height = 100 #796
terrain_width = 100 #1002
terrain = Configuration(terrain_width, terrain_height)
my_path = Path(__file__).parent / "../data/terrain/GIS_terrain_resize.csv"
terrain.load_from_csv(my_path)

# SENSOR LIST
sensor_rad = [25, 25, 15, 10]
sensor_type = ["acoustic","acoustic","acoustic","acoustic"]
num_sensors = len(sensor_type)
sensor_comm_ratio = 1.75 # ratio of sensor communication to sensing radius 
meas_type = "bearing"
# ------------------------------------------

sensor_locs = [41, 47, 49, 52]
target_loc = [45, 50]

FIM = build_FIM(sensor_locs, target_loc)
print("Baseline determinant:", np.linalg.det(FIM))


add_sensor_rad = 25
add_sensor_type = "seismic"
det = [[0 for i in range(terrain_height)] for j in range(terrain_width)]
sensor_locs.append(0) #add x
sensor_locs.append(0) #add y

for i in range(terrain_height): #height == y?
    for j in range(terrain_width): #width == x?
        # check if within range of target
        #if ( (j-target_loc[0])**2 + (i-target_loc[1])**2)**(1/2) > add_sensor_rad:
        #    det[i][j] = 0
        #else:
            sensor_locs[4], sensor_locs[5] = j, i #reset the positions
            FIM = build_FIM(sensor_locs, target_loc)
            det[i][j] = np.linalg.det(FIM)

plt.plot(sensor_locs[0], sensor_locs[1], "bx")
plt.plot(sensor_locs[2], sensor_locs[3], "bx")
plt.plot(target_loc[0], target_loc[1], "rx")


plt.imshow(det, cmap='viridis', interpolation='nearest')
plt.colorbar()
plt.show()







