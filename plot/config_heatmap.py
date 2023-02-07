# This is a python script for generating a heatmap for visualization of different configurations resulst
# The 'times' contains the binary coded uint64 vector of query times
# All input vectors are binary coded

# Example run configs:
# config_heatmap.py -c {folder_path}/configs -p {folder_path}/phases -t {folder_path}/times


import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import pandas as pd
import argparse


# command line interface
def parse_cmd():
    arg_parser = argparse.ArgumentParser()
    arg_parser.add_argument('-c', '--configs', help='file containing list of tested configs (i.e. row names)')
    arg_parser.add_argument('-p', '--phases', help='file containing list of phase names (i.e. column names')
    arg_parser.add_argument('-t', '--times', help='file containing runtime results for each configuration and phase')
    arg_parser.add_argument('-d', '--diagram', default='CRP runtimes (microsec)', help='diagram name')
    args = arg_parser.parse_args()
    print("args=%s" % args)
    return args


# helper method to read binary encoded files
def read_binary(path, type):
    if type == str:
        # reading binary string is not that easy
        with open(path, 'rb') as f:
            data = f.read()
        # split data by the null character to separate strings
        strings = [x.decode('utf-8') for x in data.split(b'\0')[:-1]]
        return strings
    else:
        with open(path, mode='rb') as file:
            numbers = np.fromfile(file, dtype=type)
        return numbers



def create_data(path_times, num_cols):
    data_1d = np.array(read_binary(path_times, np.int64))
    data_2d = data_1d.reshape(-1, num_cols)
    return data_2d


# make heatmap plot
def plot_heatmap(data, phases, configs, name):
    # set plot size
    # plt.figure(figsize=(16, 7))

    # plot the heatmap
    plt.imshow(data, cmap='hot')

    # add a color bar
    plt.colorbar()

    # set x and y labels
    plt.xlabel('Measured phases')
    plt.ylabel('Run configurations')
    plt.title(name)

    # set row, column names
    columns = phases
    rows = configs
    plt.xticks(np.arange(len(columns)), columns)
    plt.yticks(np.arange(len(rows)), rows)

    # show plot
    plt.savefig("configs_heatmap.png")
    plt.show()


# main program
args = parse_cmd()
phases = read_binary(args.phases, str)
configs = read_binary(args.configs, str)
data = create_data(args.times, len(phases))
plot_heatmap(data, phases, configs, args.diagram)
