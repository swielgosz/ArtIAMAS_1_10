from landscape import Configuration, Node
from sensors import Sensor
import numpy as np
import pprint as pp
import os
import shutil
from pathlib import Path


MAX_ITERS = 100000 # max number of iterations for Monte Carlo simulation
sensor_comm_ratio = 1.5 # ratio of sensor communication to sensing radius (because communication radius is typically larger)
CASE = 4 # case number describing combination of sensor radii and modalities:
         # !! Change later to make more general
         # Four cases for sensor radius and modality combinations:
         # Case 1: Uniform sensor modality, sensor sensing radius, and sensor communication radius
         # Case 2: Two sensor modalities (acoustic and seismic) which are uniformly distributed, and uniform sensor communication radius
         # Case 3: Uniform sensor modality, random (but bounded) sensor communication radius 
         # Case 4: Mixed sensor modalities, random sensor communication radius

# Function places sensors in unique random valid locations (i.e. not in water etc.) on map equal in size to provided terrain map
def make_basic_seismic_map(num_sensors, sensor_radius, max_width, max_height, terrain):
    sensor_map = []
    already_placed = []

    for k in range(num_sensors):

        # Place sensors in random locations that are not already occupied by another sensor
        # Terrain at location (i,j) must be suitable to be occupied by sensor (as noted by attribute is_valid_sensor_location)
        item = {}

        j = np.random.randint(0, max_width)
        i = np.random.randint(0, max_height)
        
        while (i, j) in already_placed or not terrain.is_valid_sensor_location(i, j):
            j = np.random.randint(0, max_width)
            i = np.random.randint(0, max_height)

        item["j"] = j
        item["i"] = i

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

# TASK double check what try is for
try:
    # Define number of sensors and radii
    num_sensors = 4
    if CASE == 1 or CASE == 2:
        sensor_radius = 20
    else:
        sensor_radius = np.random.uniform(low = 10, high = 20, size = num_sensors)
    for i in range(MAX_ITERS):
        print("Iteration: ", i, end="\r")
        
        # Assign sensors to grid points
        _map = make_basic_seismic_map(num_sensors, sensor_radius, terrain_width, terrain_height, terrain)

        # Calculate cost and check if configuration is valid
        score = terrain.get_configaration_observability(_map)
        valid = terrain.is_configuration_valid(_map)

        # Record each valid configuration that has a cost greater than the previous max
        if score > best_score and valid:
            print("New best score: ", score)
            best_score = score
            best_config = _map
            config_history.append(_map)
    # print(best_score)
    # pp.pprint(best_config)
            
    # Save each new best configuration
    dirpath = './case_' + str(CASE)
    if os.path.exists(dirpath) and os.path.isdir(dirpath):
        shutil.rmtree(dirpath)
        os.mkdir(dirpath)
    else:
        os.mkdir(dirpath)
    terrain.save_history(config_history,dirpath)
    terrain.plot_grid(best_config)

# Save current best configuration if user stops before max iteration
except KeyboardInterrupt:
    print("Stopped early")
    print(best_score)
    pp.pprint(best_config)
    dirpath = './case_' + str(CASE) 
    if os.path.exists(dirpath) and os.path.isdir(dirpath):
        shutil.rmtree(dirpath)
        os.mkdir(dirpath)
    else:
        os.mkdir(dirpath)
    terrain.save_history(config_history,dirpath)
    terrain.plot_grid(best_config)

