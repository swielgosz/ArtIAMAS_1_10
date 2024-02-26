from landscape import Configuration, Node
from sensors import Sensor
import numpy as np
import pprint as pp
import os
import shutil
from pathlib import Path
from scipy.optimize import minimize
from scipy import optimize
import time

sensor_comm_ratio = 1.5

# Fcn accepts a list of sensor radii, corresponding type, and location
# The sensors in the list are placed on the map
def make_basic_seismic_map(num_sensors, sensor_rad, sensor_type, x):
    sensor_map = []
    already_placed = []

    for k in range(num_sensors):

        item = {}

        # pull data from lists
        item["j"] = x[1+2*k]
        item["i"] = x[0+2*k]
        item["sensor_type"] = sensor_type[k]
        item["sensor_radius"] = sensor_rad[k]
        item["sensor_comm_radius"] = sensor_comm_ratio*sensor_rad[k]

        # Record placed sensors
        sensor_map.append(item)
    
    # Return list of sensors
    return sensor_map

terrain_height = 600
terrain_width = 600
terrain = Configuration(terrain_width, terrain_height)
my_path = Path(__file__).parent / "../data/terrain/GIS_terrain_cropped.csv"
terrain.load_from_csv(my_path)
num_sensors = 4
sensor_rad = [60,60,60,60]
# sensor_rad = [5,5,5,5]
sensor_type = ["seismic", "acoustic", "seismic", "acoustic"]
# simulated annealing
# x_final = [6*29.633102152761694,6*32.394805841102084,6*24.877332369525988,6*41.920254394788856,6*49.714748241100835,6*44.71109733315723,6*41.374726681728134,6*32.30238182521251]

x_final = [6*37,6*47,6*40,6*41,6*16,6*32,6*23,6*45]
#Monte Carlo
_map = make_basic_seismic_map(num_sensors, sensor_rad, sensor_type, x_final)
terrain.plot_grid(_map)