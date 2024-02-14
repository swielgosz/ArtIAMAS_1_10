

from landscape import Configuration, Node
from sensors import Sensor
import numpy as np
import pprint as pp
import os
import shutil
from pathlib import Path
import csv
import os.path


MAX_ITERS = 100000
sensor_comm_ratio = 1.5
CASE = 4

def make_basic_seismic_map(num_sensors, sensor_radius, max_width, max_height, terrain):
    sensor_map = []
    already_placed = []

    for k in range(num_sensors):
        item = {}

        j = np.random.randint(0, max_width)
        i = np.random.randint(0, max_height)
        
        while (i, j) in already_placed or not terrain.is_valid_sensor_location(i, j):
            j = np.random.randint(0, max_width)
            i = np.random.randint(0, max_height)

        item["j"] = j
        item["i"] = i

       
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

        
        already_placed.append((i, j))
        sensor_map.append(item)
    
    return sensor_map


terrain_width = 100
terrain_height = 100
terrain = Configuration(terrain_width, terrain_height)
my_path = Path(__file__).parent / "../data/terrain/terrain2.csv"
terrain.load_from_csv(my_path)


best_score = 0
best_config = None

config_history = []

try:
    num_sensors = 4
    if CASE == 1 or CASE == 2:
        sensor_radius = 20
    else:
        sensor_radius = np.random.uniform(low = 10, high = 20, size = num_sensors)
    for i in range(MAX_ITERS):
        print("Iteration: ", i, end="\r")

        
        _map = make_basic_seismic_map(num_sensors, sensor_radius, terrain_width, terrain_height, terrain)




        score = terrain.get_configaration_observability(_map)
        valid = terrain.is_configuration_valid(_map)

        if score > best_score and valid:
            print("New best score: ", score)
            best_score = score
            best_config = _map
            config_history.append(_map)
    # print(best_score)
    # pp.pprint(best_config)
    dirpath = './case_' + str(CASE)
    if os.path.exists(dirpath) and os.path.isdir(dirpath):
        shutil.rmtree(dirpath)
        os.mkdir(dirpath)
    else:
        os.mkdir(dirpath)
    terrain.save_history(config_history,dirpath)
    terrain.plot_grid(best_config)

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

