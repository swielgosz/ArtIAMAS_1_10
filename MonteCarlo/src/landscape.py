import matplotlib.pyplot as plt
import numpy as np
import pprint as pp


# JM - added decay to get_node_seismic_score() 
    # and get_node_acoustic_score() to Configuration class

SENSOR_RADIUS = 2 # max sensor radius
K = 10000

# For landscape classes:
# Color: used for plotting
# Land_type: type of terrain
# is_sensor_valid: describes if a sensor can be placed on that terrain

# Water class
class Water:
    def __init__(self):
        self.color = (10, 62, 235)
        self.land_type = 'water'
        self.is_sensor_valid = False 

# Land class
class Land: 
    def __init__(self):
        # self.color = (250, 245, 145)
        # self.color = (205, 235, 120)
        self.color = (235, 250, 195)
        self.land_type = 'land'
        self.is_sensor_valid = True

# Road class
class Road:
    def __init__(self):
        self.color = (0, 0, 0)
        self.land_type = 'road'
        self.is_sensor_valid = False

# Node class
class Node:
    def __init__(self, x, y, land_type=None, veg_level=0):

        self.x = x
        self.y = y

        self.land_type = land_type
        self.veg_level = veg_level

        self.sensor_type = None
        self.sensor_radius = 0


class Configuration:

    # Constructor
    def __init__(self, width, height):

        self.width = width
        self.height = height

        self.grid = []

        for i in range(height):
            row = []
            for j in range(width):
                row.append(Node(i, j, Land()))
            self.grid.append(row)

    # Return terrain type of each grid point (as a word rather than numeric value)
    def print_grid(self):
        for row in self.grid:
            for node in row:
                print(node.land_type, end=' ')
            print()

    # Load terrain file from csv
    def load_from_csv(self, filepath):
        
        # Read in numberic data (where each number in the csv file correlates to a terrain type)
        with open(filepath, 'r') as f:
            lines = f.readlines()
            for i in range(len(lines)):
                lines[i] = lines[i].strip().split(',')
                    
        # Define terrain type (and veg_level for land) according to numeric value at each grid point
        for i in range(len(lines)):
            for j in range(len(lines[i])):
                item = lines[i][j].strip()
                if item == '2':
                    self.grid[i][j].land_type = Water()
                elif item == '0':
                    self.grid[i][j].land_type = Road()
                else:
                    self.grid[i][j].land_type = Land()
                    veg_level = float(item.strip())
                    if veg_level == 1:
                        veg_level = 0
                    elif veg_level == 3:
                        veg_level = 1
                    self.grid[i][j].veg_level = veg_level

    # Check if a sensor can be placed at a specific grid point
    def is_valid_sensor_location(self, i, j):

        if self.grid[i][j].land_type.is_sensor_valid:
            return True
        
        return False
    
    # Populate grid point with sensor
    def deploy_sensor(self, i, j, sensor_type, sensor_radius):
        self.grid[i][j].sensor_type = sensor_type # "acoustic" or "seismic"
        self.grid[i][j].sensor_radius = sensor_radius

    # Clear sensors from grid points
    def clear_sensors(self):
        for i in range(self.height):
            for j in range(self.width):
                self.grid[i][j].sensor_type = None

    # PLotting the terrain with sensors
    def plot_grid(self, s_map=None):
        
        data = self.grid

        # First plot terrain with unique colors for each terrain type
        for i in range(self.height):
            for j in range(self.width):

                if self.grid[i][j].land_type.land_type == "land":
                    # Distinguish between low and high vegetation
                    # TASK change this, this was a last minute addition and not well done
                    if self.grid[i][j].veg_level < 0.5:
                        data[i][j] = self.grid[i][j].land_type.color
                    else:
                        data[i][j] = (30, 95, 5)
                elif self.grid[i][j].land_type != None:
                    data[i][j] = self.grid[i][j].land_type.color
                else:
                    data[i][j] = (255, 255, 255) # black if there is no terrain

        fig, ax = plt.subplots()

        plt.title(f"Best Configuration ")
        ax.imshow(data)
        
        # Plot sensors with unique color for each modality
        if s_map != None:
            for item in s_map:
                if item["sensor_type"] == "seismic":
                        sensor_color = '#8a54df' # purple
                elif item["sensor_type"] == "acoustic":
                    # sensor_color = '#6e5502' #'#e88e10' # orange
                    sensor_color = '#e88e10' # orange

                # Plotting sensors location
                circle = plt.Circle((item["j"], item["i"]), item["sensor_radius"], color = sensor_color, linewidth = 2, fill=False)
                ax.add_patch(circle)

                # Plotting sensor sensing radius
                circle2 = plt.Circle((item["j"], item["i"]), .5, color = sensor_color,  fill=True)
                ax.add_patch(circle2)

                # Plotting sensor communication radius
                circle3 = plt.Circle((item["j"], item["i"]), item["sensor_comm_radius"], color = sensor_color,linestyle='--', linewidth = 0.5,fill=False)
                ax.add_patch(circle3)

                

        plt.show()
        

    # Save plots. This is basically the same as plot_grid. Not efficient
    def save_history(self, history,dirpath):
        c = 0
        for s_map in history:
            data = [row[:] for row in self.grid]

            for i in range(self.height):
                for j in range(self.width):

                    if self.grid[i][j].land_type.land_type == "land":
                        if self.grid[i][j].veg_level < 0.5:
                            data[i][j] = self.grid[i][j].land_type.color
                        else:
                            data[i][j] = (30, 95, 5)
                    elif self.grid[i][j].land_type != None:
                        data[i][j] = self.grid[i][j].land_type.color
                    else:
                        data[i][j] = (255, 255, 255)

            fig, ax = plt.subplots()

            
            ax.imshow(data)
            
            if s_map != None:
                for item in s_map:
                    if item["sensor_type"] == "seismic":
                        sensor_color = '#8a54df' # purple
                    elif item["sensor_type"] == "acoustic":
                        # sensor_color = '#6e5502' #
                        sensor_color = '#e88e10' # orange

                    circle = plt.Circle((item["j"], item["i"]), item["sensor_radius"], color = sensor_color, linewidth = 2, fill=False)
                    ax.add_patch(circle)

                    circle2 = plt.Circle((item["j"], item["i"]), .5, color = sensor_color, fill=True)
                    ax.add_patch(circle2)

                    circle3 = plt.Circle((item["j"], item["i"]), item["sensor_comm_radius"], color = sensor_color, linestyle='--', linewidth = 0.5, fill=False)
                    ax.add_patch(circle3)

                
            plt.title(f"Best Configuration {c+1}")
            plt.savefig(f"{dirpath}/best_config_{c}.png")

            c += 1


    # Assign acoustic sensor component of "coverage score" (defined in project PPT)
    # JM - add "gaussian" ability here
        # pass in coords, distance from sensor, and std dev
        # Upping sigma would increase the "decay" of the sensing ability    
    def get_node_acoustic_score(self, i, j, d, sigma):
        
        n = self.grid[i][j]
        a_0 = 1

        if n.land_type.land_type == "water":
            scale_factor = 0.5
        elif n.land_type.land_type == "land":
            scale_factor = 1 / (n.veg_level + 1) # accounts for low vs high vegetation
        elif n.land_type.land_type == "road":
            scale_factor = 1
        G = np.exp(-0.5*(d/sigma)**2)
        # return something like a_0 * scale * gaussian scaling (exists in [0, 1])
        return a_0 * scale_factor * G

    # Assign seismic sensor component of "coverage score" (defined in project PP
    def get_node_seismic_score(self, i, j, d, sigma):
        
        n = self.grid[i][j]
        a_0 = 1

        if n.land_type.land_type == "water":
            scale_factor = .25
        elif n.land_type.land_type == "land":
            scale_factor = 1
        elif n.land_type.land_type == "road":
            scale_factor = 1

        G = np.exp(-0.5*(d/sigma)**2)
        # return something like a_0 * scale * gaussian scaling (exists in [0, 1])
        return a_0 * scale_factor * G
    
    # Weighting factor of each terrain type (alpha in PPT)
    def get_desirability_score(self, i, j):
        
        n = self.grid[i][j]

        if n.land_type.land_type == "water":
            return .2
        elif n.land_type.land_type == "land":
            return 1
        elif n.land_type.land_type == "road":
            return 2
    
    # Calculate cost (configuration profit in PPT) of configuration. We want to maximize the cost, i.e. maximize information known about the environment
    def get_configaration_observability(self, s_map=None):
        if s_map == None:
            s_map = self.get_sensor_map()
        
        conf_score = 0

        # Loop through each grid point
        for i in range(self.height):
            for j in range(self.width):
                
                n_acoustic = 0
                n_seismic = 0

                # Check which sensors each grid point is detected by, and assign "coverage score" accordingly
                for sensor in s_map:
                    # JM - will reuse the distance calculation for gaussian decay
                    dist = np.sqrt((i - sensor["i"])**2 + (j - sensor["j"])**2)
                    if  dist <= sensor["sensor_radius"]:

                        if sensor["sensor_type"] == "acoustic" and n_acoustic == 0:
                            n_acoustic = self.get_node_acoustic_score(i, j, dist, sigma=sensor['sensor_radius'])

                        elif sensor["sensor_type"] == "seismic" and n_seismic == 0:
                            n_seismic = self.get_node_seismic_score(i, j, dist, sigma=sensor['sensor_radius'])

                # Calculate "profit" at each grid point by weighting "coverage score" by terrain weighting factor, and sum to find total "configuration profit"
                conf_score += (n_acoustic + n_seismic) * self.get_desirability_score(i, j)

        return conf_score
    
    # Return informaton about sensors in map
    def get_sensor_map(self):

        sensor_map = []
        sensor_count = 0

        for i in range(self.height):
            for j in range(self.width):
                if self.grid[i][j].sensor_type != None:
                    sensor_map.append({
                        "i": i,
                        "j": j,
                        "sensor_type": self.grid[i][j].sensor_type,
                        "sensor_radius": self.grid[i][j].sensor_radius,
                        "sensor_id": sensor_count
                    })
                sensor_count += 1
                    
        
        return sensor_map

    # Check if the sensors in a configuration are connected via a WSN
    def is_configuration_valid(self, s_map=None):
        if s_map == None:
            s_map = self.get_sensor_map()

        graph = {}

        # For each sensor, record which sensors are within communication radius
        for i in range(len(s_map)):
            curr_node_connections = []
            for j in range(len(s_map)):
                if i == j:
                    continue
                else:
                    if np.sqrt((s_map[i]["i"] - s_map[j]["i"])**2 + (s_map[i]["j"] - s_map[j]["j"])**2) <= s_map[j]["sensor_comm_radius"]:
                        curr_node_connections.append(str(j))
            graph[str(i)] = curr_node_connections
                

        # Implement breadth first search to check if sensors are connected
        visited = [] # List for visited nodes.
        queue = []     #Initialize a queue
        visited.append('0')
        queue.append('0')

        while queue:          # Creating loop to visit each node
            m = queue.pop(0) 

            for neighbour in graph[m]:
                if neighbour not in visited:
                    visited.append(neighbour)
                    queue.append(neighbour)
        
        # If list of visited neighbors is equal in length to number of sensors, sensors are connected
        if len(visited) == len(s_map):
            return True

if __name__ == '__main__':
    landscape = Configuration(100, 100)