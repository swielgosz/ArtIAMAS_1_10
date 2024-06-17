def make_basic_seismic_map(num_sensors, sensor_rad, sensor_type, sensor_mode, sensor_comm_ratio, x):
    sensor_map = []
    for k in range(len(x)//2):

        item = {}

        # pull data from lists
        item["j"] = x[1+2*k]
        item["i"] = x[0+2*k]
        item["sensor_type"] = sensor_type[k]
        item["sensor_radius"] = sensor_rad[k]
        item["sensor_mode"] = sensor_mode[k]
        item["sensor_comm_radius"] = sensor_comm_ratio*sensor_rad[k]

        # Record placed sensors
        sensor_map.append(item)
    
    # Return list of sensors
    return sensor_map