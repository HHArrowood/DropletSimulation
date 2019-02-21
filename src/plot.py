import csv
import numpy as np
import matplotlib.pyplot as plt
import os

def get_interface_data(filename):
    with open(filename, 'r') as f:
        reader = csv.reader(f)
        in_header = True
        x = []
        y = []
        for row in reader:
            if in_header:
                in_header = False
                continue
            x.append(float(row[0]))
            y.append(float(row[1]))
    return x, y
            
def plot_interface(filename, plot_filename):
    x, y = get_interface_data(filename)
    plt.figure()
    plt.plot(x,y,'.')
    plt.xlabel("x")
    plt.ylabel("y")
    plt.savefig(plot_filename)
    plt.close()

fnames = [f for f in os.listdir('.') if os.path.isfile(f)]
for fname in fnames:
    if fname.startswith("interface") and fname.endswith(".dat"):
        i_dat = fname.index(".dat")
        raw_name = fname[:i_dat]
        plot_fname = raw_name + ".png"
        plot_interface(fname, plot_fname)
