from .targets import target

def add_targets(_sensor_map, target_locs):
    sensor_map = []
    for k in range(len(target_locs)//2):

        _sensor_map.plot(target_locs[0+2*k], target_locs[1+2*k], "r^")
        # _sensor_map.scatter(target_locs[0+2*k], target_locs[1+2*k], c = target.color)
    
    # Return the figure
    return _sensor_map


def add_targets_localized(_sensor_map, target_locs, localized):
    sensor_map = []
    for k in range(len(target_locs)//2):
        if localized[k] == 1:
            _sensor_map.plot(target_locs[0+2*k], target_locs[1+2*k], "k^")
        if localized[k] == 0:
            _sensor_map.plot(target_locs[0+2*k], target_locs[1+2*k], "r^")
        # _sensor_map.scatter(target_locs[0+2*k], target_locs[1+2*k], c = target.color)
    
    # Return the figure
    return _sensor_map