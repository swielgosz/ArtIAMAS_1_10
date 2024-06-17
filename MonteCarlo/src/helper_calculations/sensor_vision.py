import numpy as np
from .landscape import Configuration, Node

def sensor_target_angle(sensor, target):
    dx, dy = (target[0] - sensor[0]), (target[1] - sensor[1])
    return np.arctan(dy/dx)

def get_LOS_coeff(sensor, target, terrain):
    '''
    Returns the LOS coefficient that scales the sensor noise

    INPUTS: One sensor location, one target and the map
    RETURNS: scaling value
    '''
    dx, dy = (target[0] - sensor[0]), (target[1] - sensor[1])
    unit_steps = int((dx**2+dy**2)**(1/2))
    if np.isclose(dx, 0, rtol=1e-06) == True:
        theta = np.pi
    else: 
        theta = np.arctan2(dy, dx)
    
    LOS_coeff = 0
    for k in range(unit_steps):
        i = round(sensor[0] + (k+1)*np.cos(theta))
        j = round(sensor[1] + (k+1)*np.sin(theta))
        LOS_coeff+=terrain.LOS_loss_coeff(j, i)

    if unit_steps != 0:
        return LOS_coeff/unit_steps
    else:
        return 1

def sensor_vision(sensor, sensor_rad, target, meas_type):
    '''
    Returns a set of nodes that the sensor sees for plotting purposes
    EX: Bearing will return a straight line that gives the potential locations of the target
    EX: Distance will return a circle of potential locs

    INPUTS: One sensor location, its radius, the target loc, and the type
    RETURNS: points to plot that correspond to the sensor type
    '''
    
    nodes = []
    # Flag if the target doesn't see the sensor
    detected = False

    # first check if the target is actually in the sensor radius
    dist = ((sensor[0] - target[0])**2 + (sensor[1] - target[1])**2)**(1/2)
    #print(dist, sensor, sensor_rad)
    if dist <= sensor_rad:
        detected = True
        if meas_type == "bearing":
            # Finds the bearing from the sensor to target
            dx, dy = (target[0] - sensor[0]), (target[1] - sensor[1])
            if np.isclose(dx, 0, rtol=1e-06) == True:
                theta = np.pi
            else: 
                theta = np.arctan2(dy, dx)

            # Compensate for negative/negative calcelation
            #if dx < 0 and dy < 0:
            #    theta += np.pi
            
            # Calculates the nodes seen in this line of sight
            unit_vec = [np.cos(theta), np.sin(theta)]
            nodes = [alpha*np.array(unit_vec) + np.array(sensor) for alpha in range(1,sensor_rad)]
            
            # Puts it into a list and converts to integers for plotting
            nodes = [node.tolist() for node in nodes]
            #nodes = [[int(cord) for cord in node] for node in nodes]
            #print(nodes)

        if meas_type == "radial":
            nodes = []
            # Given a calculated distance, return the nodes encircling
            n_pts = 50
            nodes_x = [dist*np.sin(th)+sensor[0] for th in np.linspace(0, 2*3.14, n_pts)]
            nodes_y = [dist*np.cos(th)+sensor[1] for th in np.linspace(0, 2*3.14, n_pts)]

            nodes = [(nodes_x[i], nodes_y[i]) for i in range(0, n_pts)]
            # TO DO: Calcs here for sensor vision

        return nodes, detected

    # return an empty list otherwise
    else:
        return nodes, detected

def sensor_reading(sensor, sensor_rad, target, meas_type):
    '''
    Returns the measurement from the sensors
    Bearing: will return an (x, y) unit vector in direction from sensor to target
    Radii: will return a distance

    INPUTS: One sensor location, its radius, the target loc, and the type
    RETURNS: Either a unit vector tuple or a single value
    '''

    # first check if the target is actually in the sensor
    dist = ((sensor[0] - target[0])**2 + (sensor[1] - target[1])**2)**(1/2)

    detected = False
    if dist <= sensor_rad:
        detected = True
        # BEARING: return the unit vec
        if meas_type == "bearing":
            return_meas = []
            # Finds the bearing from the sensor to target
            dx, dy = (target[0] - sensor[0]), (target[1] - sensor[1])
            theta = np.arctan(dy/dx)
            
            # Compensate for negative/negative calcelation
            if dx < 0 and dy < 0:
                theta += np.pi
            
            # Calculates the nodes seen in this line of sight
            return_meas = [np.cos(theta), np.sin(theta)]
    
        # DISTANCE: Return the radius
        if meas_type == "radial":
            return_meas = dist
    

    return (return_meas, detected)

def target_localization_bearing(bearing_units, sensors):
    '''
    Localizes a target with two bearing sensors 

    INPUTS: The bearing unit vectors (list of lists), sensor locations (list)
    RETURNS: The target location as a list
    '''
    # For 2 bearing sensors, return the intersecting point
    m = [] #set empty vector for slope calcs
    
    # Construct the list of slopes
    for i, bearing in enumerate(bearing_units):
        dy, dx = bearing[1], bearing[0]
        m.append(dy/dx)

    # Now construct the x, y list
    x, y = [], []
    for i in range(len(sensors)//2):
        xi, yi = sensors[0+2*i], sensors[1+2*i]
        x.append(xi)
        y.append(yi)

    print(x,y)

    # Perform the intersection calculation
    x_star = (m[0]*x[0] - m[1]*x[1] + y[1] - y[0])/(m[0]-m[1])
    y_star = m[0]*(x_star-x[0]) + y[0]
        
    return [x_star, y_star]

def target_localization_radius(radii, sensors):
    '''
    Localizes a target with two bearing sensors

    INPUTS: Three radii (list), sensor positions (list)
    RETURNS: The target location as a list
    '''

    # The sensor only detects a radius
    # Therefore, two sensors can localize a target to two points

    # Start by finding the intersection of two circles
    # from https://math.stackexchange.com/questions/256100/how-can-i-find-the-points-at-which-two-circles-intersect
    
    x, y, r = [], [], []
    for i in range(2):
        xi, yi = sensors[0+2*i], sensors[1+2*i]
        x.append(xi)
        y.append(yi)
        r.append(radii[i])
        
    # Perform the calculations
    R = ((x[1]-x[0])**2+(y[1]-y[0])**2)**(1/2)
    sqrt_term = 0.5*( 2*(r[0]**2+r[1]**2)/(R**2) - (r[0]**2-r[1]**2)**2/(R**4) - 1)**(1/2)

    x_inter1 = 0.5*(x[0]+x[1])+((r[0]**2-r[1]**2)/(2*R**2))*(x[1]-x[0])+sqrt_term*(y[1]-y[0])
    y_inter1 = 0.5*(y[0]+y[1])+((r[0]**2-r[1]**2)/(2*R**2))*(y[1]-y[0])+sqrt_term*(x[0]-x[1])

    x_inter2 = 0.5*(x[0]+x[1])+((r[0]**2-r[1]**2)/(2*R**2))*(x[1]-x[0])-sqrt_term*(y[1]-y[0])
    y_inter2 = 0.5*(y[0]+y[1])+((r[0]**2-r[1]**2)/(2*R**2))*(y[1]-y[0])-sqrt_term*(x[0]-x[1])

    # To get the actual sensor position, we need three circles to intersect
    # Using the third sensor, check which point satisfies the third condition
        # (x_sensor3 - x_inter_i)^2+(y_sensor3 - y_inter_i)^2 == r3^2

    x3, y3 = sensors[4], sensors[5]

    d1 = ((x3 - x_inter1)**2 + (y3 - y_inter1)**2)**(1/2)
    d2 = ((x3 - x_inter2)**2 + (y3 - y_inter2)**2)**(1/2)

    if abs(d1-radii[2]) < abs(d2-radii[2]):
        # If satisfies candidate 1 conditions:
        return [x_inter1, y_inter1]
    else: # if satisfies candidate 2
        return [x_inter2, y_inter2]
    

def target_localization_both(measurements, sensors, sensor_types):
    '''
    Localizes a target with two bearing sensors

    INPUTS: Radii sensor radius, bearing sensor bearing, and list of positions
    RETURNS: The target location as a list
    '''
    # Calculated from: https://math.stackexchange.com/questions/228841/how-do-i-calculate-the-intersections-of-a-straight-line-and-a-circle
    sensor_b, sensor_s =[], []
    for i in range(len(sensors)//2):
        if sensor_types[i] == "bearing":
            sensor_b.append(sensors[0+2*i])
            sensor_b.append(sensors[1+2*i])
            bearing = measurements[i]
        if sensor_types[i] == "radial":
            sensor_s.append(sensors[0+2*i])
            sensor_s.append(sensors[1+2*i])
            r = measurements[i]

    dy, dx = bearing[1], bearing[0]
    m = dy/dx
    c = -m*sensor_b[0]+sensor_b[1]
    ys, xs = sensor_s[1], sensor_s[0]

    A = m**2+1
    B = 2*(m*c - m*ys - xs)
    C = (ys**2 - r**2 + xs**2 - 2*c*ys + c**2)

    x_int1, x_int2 = (-B+(B**2-4*A*C)**(1/2))/(2*A), (-B-(B**2-4*A*C)**(1/2))/(2*A)
    y_int1, y_int2 = m*x_int1+c, m*x_int2+c

    # Pick the one that satisfies the sensing radius condition
    # on the bearing sensor (i.e. the lesser distance one)

    dist1 = ((x_int1-sensor_b[0])**2 + ((y_int1-sensor_b[1])**2) )**(1/2)
    dist2 = ((x_int2-sensor_b[0])**2 + ((y_int2-sensor_b[1])**2) )**(1/2)

    if dist1 < dist2:
        return [x_int1, y_int1]
    else:
        return [x_int2, y_int2]