# This script saves off data that we can plot later
import os

def save_data(data, names, file_name, description):
    '''
    This function will save off data in the following format as a csv:
    -------
    File description
    var 1 # # # # # ...
    var 2 # # # ....
    ...
    -------
    ACCEPTS: data - list of lists for the variable series saved off
             names - associated names for the variables
             file_name - desired name of the csv
             description - header
    RETURNS: nothing, saves data in the saved_data folder
    '''

    # Check that all saved variables have a name
    if len(data) != len(names):
        print("Need an equal set of names and data to save!")
        print("No data saved. Continuing simulation")
    
    
    # If they do, continue
    else:
        # Generate file name
        full_file_name = "./Optimization_alg/src/Saved_data/" + file_name + ".csv"
        
        # Save off data
        with open(full_file_name, "w+") as fi:
            
            # Add header identifying the file
            fi.write(description)
            fi.write("\n")
            
            # Loop over data
            for j, data_series in enumerate(data):
                fi.write(f"{names[j]}: ") # Write the name first
                
                # Now write the series
                for i in range(len(data_series)):
                    fi.write(f"{data_series[i]} ")
                
                # Move to new line
                fi.write(f"\n")
