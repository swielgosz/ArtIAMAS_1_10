from MonteCarlo.src.landscape import Configuration, Node
from sensors import Sensor
import numpy as np
import pprint as pp
import os
import shutil
from pathlib import Path
from scipy.optimize import minimize
import matplotlib.pyplot as plt
import seaborn as sns


sensor_comm_ratio = 1.5 # ratio of sensor communication to sensing radius (because communication radius is typically larger)
CASE = 1 # case number describing combination of sensor radii and modalities:


# Function places sensors in unique random valid locations (i.e. not in water etc.) on map equal in size to provided terrain map
def make_basic_seismic_map(num_sensors, sensor_radius, max_width, max_height, terrain, x):
    sensor_map = []
    already_placed = []

    for k in range(num_sensors):

        # Place sensors in random locations that are not already occupied by another sensor
        # Terrain at location (i,j) must be suitable to be occupied by sensor (as noted by attribute is_valid_sensor_location)
        item = {}

        item["j"] = x[1]
        item["i"] = x[0]

        # !! Change later to make more general
        # Assign sensor radius and modality according to case 
        if CASE == 1:
            item["sensor_type"] = "seismic"
            item["sensor_radius"] = sensor_radius
            item["sensor_comm_radius"] = sensor_comm_ratio*sensor_radius
        elif CASE == 2:
            if k <= num_sensors/2 - 1:
                item["sensor_type"] = "seismic"
            else:
                item["sensor_type"] = "acoustic"
            item["sensor_radius"] = sensor_radius
            item["sensor_comm_radius"] = sensor_comm_ratio*sensor_radius
        elif CASE == 3:
            item["sensor_type"] = "seismic"
            item["sensor_radius"] = sensor_radius[k]
            item["sensor_comm_radius"] = sensor_comm_ratio*sensor_radius[k]
        elif CASE == 4:
            if k <= num_sensors/2 - 1:
                item["sensor_type"] = "seismic"
            else:
                item["sensor_type"] = "acoustic"
            item["sensor_radius"] = sensor_radius[k]
            item["sensor_comm_radius"] = sensor_comm_ratio*sensor_radius[k]

        # Record placed sensors
        already_placed.append((i, j))
        sensor_map.append(item)
    
    # Return list of sensors
    return sensor_map



# Define terrain 
terrain_width = 100
terrain_height = 100
terrain = Configuration(terrain_width, terrain_height)
my_path = Path(__file__).parent / "../data/terrain/terrain2.csv"
terrain.load_from_csv(my_path)

# Initialize optimal sensor configuration cost and map
best_score = 0
best_config = None

config_history = []
score = [[0 for i in range(terrain_height)] for j in range(terrain_width)]
loc = []
x = [0,0]

    # Define number of sensors and radii
num_sensors = 1
if CASE == 1 or CASE == 2:
    sensor_radius = 20


for i in range(terrain_height):
    for j in range(terrain_width):
        x = [i, j]
        loc.append(x)
        if not terrain.is_valid_sensor_location(i, j):
            score_iter = 0
            #print(score_iter)
            score[i][j] = score_iter
        else:
            _map = make_basic_seismic_map(num_sensors, sensor_radius, terrain_width, terrain_height, terrain, x)
            score_iter = terrain.get_configaration_observability(_map)
            # Calculate cost and check if configuration is valid
            #print(score_iter)
            score[i][j] = score_iter



plt.imshow(score, cmap='viridis', interpolation='nearest')
plt.show()




