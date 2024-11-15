import csv

def get_csv_size(my_path):
    with open(my_path) as f:
        data = list(csv.reader(f))
    rows = len(data)
    cols = len(data[0])
    return rows,cols