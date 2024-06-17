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

# Classes:
# 1: water
# 2: road
#3: low vegetation
# 4: High vegetation
# 5: Man-made objects

# Water class
class Water:
    def __init__(self):
        # self.color = (10, 62, 235)
        self.color = (0, 0, 232)
        self.land_type = 'water'
        self.is_sensor_valid = False 

# Land class
class Land: 
    def __init__(self):
        # self.color = (250, 245, 145)
        # self.color = (205, 235, 120)
        # self.color = (235, 250, 195)
        self.color =(186, 239, 152)
        self.land_type = 'land'
        self.is_sensor_valid = True

# Road class
class Road:
    def __init__(self):
        # self.color = (0, 0, 0)
        self.color = (128, 64, 128)
        self.land_type = 'road'
        self.is_sensor_valid = False

# Manmade class
class Manmade:
    def __init__(self):
        # self.color = (168, 76, 50)
        # self.color = (239, 37, 7)
        self.color = (168, 19, 10)
        self.land_type = 'manmade'
        self.is_sensor_valid = False

# Node class
class Node:
    def __init__(self, x, y, land_type=None, veg_level=0):

        self.x = x
        self.y = y

        self.land_type = land_type
        self.veg_level = veg_level

        self.sensor_type = None # acoustic or seismic
        self.sensor_mode = None # radius or bearing
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
                # print(i,j)
            self.grid.append(row)
        # print(i,j)

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
                if item == '1':
                    self.grid[i][j].land_type = Water()
                elif item == '2':
                    self.grid[i][j].land_type = Road()
                elif item == '5':
                    self.grid[i][j].land_type = Manmade()
                else:
                    # print("here")
                    self.grid[i][j].land_type = Land()
                    veg_level = float(item.strip())
                    if veg_level == 3:
                        veg_level = 0
                    elif veg_level == 4:
                        veg_level = 1
                    self.grid[i][j].veg_level = veg_level

    # Check if a sensor can be placed at a specific grid point
    def is_valid_sensor_location(self, i, j):

        if self.grid[i][j].land_type.is_sensor_valid:
            return True
        
        return False
    
    # Populate grid point with sensor
    def deploy_sensor(self, i, j, sensor_type, sensor_mode, sensor_radius):
        self.grid[i][j].sensor_mode = sensor_mode # "bearing" or "radius"
        self.grid[i][j].sensor_type = sensor_type # "acoustic" or "seismic"
        self.grid[i][j].sensor_radius = sensor_radius

    # Clear sensors from grid points
    def clear_sensors(self):
        for i in range(self.height):
            for j in range(self.width):
                self.grid[i][j].sensor_type = None
                self.grid[i][j].sensor_mode = None

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
                        # data[i][j] = (30, 95, 5)
                        data[i][j] = (5, 129, 35)
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
                        # sensor_color = '#8a54df' # purple
                        sensor_color = '#000000' # black

                elif item["sensor_type"] == "acoustic":
                    # sensor_color = '#6e5502' #'#e88e10' # orange
                    sensor_color = '#e88e10' # orange
                    # sensor_color = '#B3ADAD' #gray
                    # sensor_color = '#FBBD1F' # more orange


                # Plotting sensors sensing radius
                circle = plt.Circle((item["i"], item["j"]), item["sensor_radius"], color = sensor_color, linewidth = 2, fill=False)
                ax.add_patch(circle)

                # Plotting sensors location
                circle2 = plt.Circle((item["i"], item["j"]), 0.5, color = sensor_color,  fill=True)
                ax.add_patch(circle2)

                # Plotting sensor communication radius
                #circle3 = plt.Circle((item["i"], item["j"]), item["sensor_comm_radius"], color = sensor_color,linestyle='--', linewidth = 0.5,fill=False)
                #ax.add_patch(circle3)
        
        return ax

        # plt.show()    

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

                    #circle3 = plt.Circle((item["j"], item["i"]), item["sensor_comm_radius"], color = sensor_color, linestyle='--', linewidth = 0.5, fill=False)
                    #ax.add_patch(circle3)

                
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
            scale_factor = 5
        elif n.land_type.land_type == "manmade":
            scale_factor = 1
        else:
            scale_factor = 0
        # G = np.exp(-0.5*(d/sigma)**2)
        # return something like a_0 * scale * gaussian scaling (exists in [0, 1])
        # return a_0 * scale_factor * G
        return a_0 * scale_factor

    # Assign seismic sensor component of "coverage score" (defined in project PPT)
    def get_node_seismic_score(self, i, j, d, sigma):
        
        n = self.grid[i][j]
        a_0 = 1

        if n.land_type.land_type == "water":
            scale_factor = .25
        elif n.land_type.land_type == "land":
            scale_factor = 1
        elif n.land_type.land_type == "road":
            scale_factor = 5
        elif n.land_type.land_type == "manmade":
            scale_factor = 0.75
        else:
            scale_factor = 0

        # G = np.exp(-0.5*(d/sigma)**2)
        # return something like a_0 * scale * gaussian scaling (exists in [0, 1])
        # return a_0 * scale_factor * G
        return a_0 * scale_factor
    
    # Weighting factor of each terrain type (alpha in PPT)
    def get_desirability_score(self, i, j):
        
        n = self.grid[i][j]

        if n.land_type.land_type == "water":
            return .1
        elif n.land_type.land_type == "land":
            return 0.5
        elif n.land_type.land_type == "road":
            return 3
        elif n.land_type.land_type == "manmade":
            return 1.5
        
    def LOS_loss_coeff(self, i, j):
        
        n = self.grid[i][j]

        if n.land_type.land_type == "water":
            return 2.5
        elif n.land_type.land_type == "land":
            if n.veg_level < 0.5:
                return 1
            elif n.veg_level > 0.5:
                return 1.5
        elif n.land_type.land_type == "road":
            return 2
        elif n.land_type.land_type == "manmade":
            return 2.5
        
    # Commented this out because we want to remove the sensor weighting factor
    # # Calculate cost (configuration profit in PPT) of configuration. We want to maximize the cost, i.e. maximize information known about the environment
    # def get_configuration_observability(self, s_map=None):
    #     if s_map == None:
    #         s_map = self.get_sensor_map()
        
    #     conf_score = 0

    #     # Loop through each grid point
    #     for i in range(self.height):
    #         for j in range(self.width):
                
    #             n_acoustic = 0
    #             n_seismic = 0

    #             # Check which sensors each grid point is detected by, and assign "coverage score" accordingly
    #             for sensor in s_map:
    #                 # JM - will reuse the distance calculation for gaussian decay
    #                 dist = np.sqrt((i - sensor["i"])**2 + (j - sensor["j"])**2)
    #                 if  dist <= sensor["sensor_radius"]:

    #                     if sensor["sensor_type"] == "acoustic" and n_acoustic == 0:
    #                         n_acoustic = self.get_node_acoustic_score(i, j, dist, sigma=sensor['sensor_radius'])

    #                     elif sensor["sensor_type"] == "seismic" and n_seismic == 0:
    #                         n_seismic = self.get_node_seismic_score(i, j, dist, sigma=sensor['sensor_radius'])

    #             # Calculate "profit" at each grid point by weighting "coverage score" by terrain weighting factor, and sum to find total "configuration profit"
    #             conf_score += (n_acoustic + n_seismic) * self.get_desirability_score(i, j) # in line below, remove n_acoustic and n_seismic so we only check if a location is covered and scale it by terrain desirability

    #     return conf_score
    
    # ADD IN LOS CALCS
    # Calculate cost (configuration profit in PPT) of configuration. We want to maximize the cost, i.e. maximize information known about the environment
    def get_configuration_observability(self, s_map=None):
        if s_map == None:
            s_map = self.get_sensor_map()
        
        conf_score = 0

        # Loop through each grid point
        for i in range(self.height):
            for j in range(self.width):

                # Check which sensors each grid point is detected by, and assign "coverage score" accordingly
                for sensor in s_map:
                    # JM - will reuse the distance calculation for gaussian decay
                    dist = np.sqrt((i - sensor["i"])**2 + (j - sensor["j"])**2)
                    if  dist <= sensor["sensor_radius"]:

                        # If a sensor sees a node, break and move on to next node
                        conf_score += self.get_desirability_score(i, j)
                        break

        return conf_score
    
    # ADD IN CALCULATION FOR LOCALIZATION SCORE
    def get_localization_score(self, s_map=None):
        if s_map == None:
            s_map = self.get_sensor_map()
        
        localization_score = 0
        
        # Loop through each grid point
        for i in range(self.height):
            for j in range(self.width):

                localized_flag = 0
                radius_sensor_count = 0
                bearing_sensor_count = 0

                bearing_sublist = []
                radius_sublist = []

                # If grid point is sensed by 3 radius sensors, 
                    # 2 bearing sensors, or 
                    # 1 radius and 1 bearing 
                    # we may be able to localize it!
                for sensor in s_map:

                    dist = np.sqrt((i - sensor["i"])**2 + (j - sensor["j"])**2)
                    if  dist <= sensor["sensor_radius"]:
                        if sensor["sensor_mode"] == "bearing":
                            bearing_sensor_count += 1
                            bearing_sublist.append([sensor["i"],sensor["j"]])

                        elif sensor["sensor_mode"] == "radius":
                            radius_sensor_count += 1
                            radius_sublist.append([sensor["i"],sensor["j"]])


                        # Count if the num of required sensors are there
                        CND1, CND2 = bearing_sensor_count >= 2, radius_sensor_count >= 3, 
                        CND3 = (bearing_sensor_count >= 1 & radius_sensor_count >= 1)

                        if CND1: # if the bearing sublist is > 2
                            # Create a list of duplicated. 
                            dups = {tuple(x) for x in bearing_sublist if bearing_sublist.count(x)>1}
                            # If empty, return TRUE
                            sub_CND1 = (dups == set())
                        
                        if CND2:# if the radius sublist is > 3
                            # Create a list of duplicated. 
                            dups = {tuple(x) for x in radius_sublist if radius_sublist.count(x)>1}
                            sub_CND2 = (dups == set())

                        # Now check that the required num of sensors meet the config reqs
                        CND1_tot = CND1 and sub_CND1
                        CND2_tot = CND2 and sub_CND2

                        # If true, this pt is satsified and break
                        if (CND1_tot|CND2_tot|CND3) == True:
                            localized_flag = 1
                            break
                
                # Compute the localization score using the {0,1} binary
                localization_score += localized_flag*self.get_desirability_score(i, j)
                        
        return localization_score
    
    # ADD IN CALCULATION FOR TOTAL MAP SCORE

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
                        "sensor_mode": self.grid[i][j].sensor_mode,
                        "sensor_radius": self.grid[i][j].sensor_radius,
                        "sensor_id": sensor_count
                    })
                sensor_count += 1
                    
        
        return sensor_map

    # Check if the sensors in a configuration are connected via a WSN
    # UPDATE THIS for 2 way communication
    def is_configuration_valid(self, s_map=None):
        #Inits
        WSN = False
        
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
                    if (np.sqrt((s_map[i]["i"] - s_map[j]["i"])**2 + (s_map[i]["j"] - s_map[j]["j"])**2) <= s_map[j]["sensor_comm_radius"]) & (np.sqrt((s_map[i]["i"] - s_map[j]["i"])**2 + (s_map[i]["j"] - s_map[j]["j"])**2) <= s_map[i]["sensor_comm_radius"]):
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
            WSN = True

        return WSN

if __name__ == '__main__':
    landscape = Configuration(100, 100)
    # landscape = Configuration(1002, 796)