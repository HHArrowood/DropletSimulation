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
    
def generate_droplet_surface():    
    t = np.linspace(-3.14/2,3.14/2,100) 
    x1 = (.05075-2*.00093)*np.cos(t)
    y1 = (.05075-2*.00093)*np.sin(t)
    return x1, y1
         
def plot_interface(filename, plot_filename):
    x, y = get_interface_data(filename)
    plt.figure()
    plt.plot(x,y,linewidth=.3)
    x1, y1 = generate_droplet_surface()
    plt.plot(x1,y1,linewidth=.3)
    plt.xlabel("x")
    plt.ylabel("y")
    plt.savefig(plot_filename, dpi = 400)
    plt.close()
    plt.axes().set_aspect('equal')

def get_slice_after_last_space(string):
    strrev = string[::-1]
    last_space_idx = strrev.index(" ")
    strrev_sliced = strrev[:last_space_idx]
    string_sliced = strrev_sliced[::-1]
    return string_sliced

fnames = [f for f in os.listdir('.') if os.path.isfile(f)]
for fname in fnames:
    if fname.startswith("interface") and fname.endswith(".dat"):
        i_dat = fname.index(".dat")
        raw_name = fname[:i_dat]
        name_no_spaces = get_slice_after_last_space(raw_name)
        plot_fname = name_no_spaces + ".png"
        plot_interface(fname, plot_fname)
