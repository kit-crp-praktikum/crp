# This is a python script for generating a speedup plot
# The 'input' contains the binary coded uint64 vector of query times
# All input vectors are binary coded

# Example run configs:
# speedup.py -c {folder_path}/configs -t {folder_path}/threads -i {folder_path}/times


import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import pandas as pd
import argparse


# command line interface
def parse_cmd():
    arg_parser = argparse.ArgumentParser()
    arg_parser = argparse.ArgumentParser()
    arg_parser.add_argument('-c', '--configs', help='file containing list of tested configs (i.e. row names)')
    arg_parser.add_argument('-t', '--threads', help='tested thread numbers for each algorithm (i.e. column names')
    arg_parser.add_argument('-i', '--input', help='file containing runtime results for each configuration and thread number')
    arg_parser.add_argument('-d', '--diagram', default='CRP speedups (microsec)', help='diagram name')

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


# create whole dataframe
def create_data(path_times, num_threads):
    times = read_binary(path_times, np.int64)
    data = []
    for i in range(0, len(times), num_threads):
        # times[i]: single-thread runtime for current config
        config_times = [times[i]/x for x in times[i:i + num_threads]]
        data.append(config_times)
    return data


# make speedup plot
# data contains speedups for each config
def plot_speedup(data, threads, configs, name):
    x = threads
    for i in range (len(data)):
        plt.plot(x, data[i], label=configs[i], marker='o')

    plt.title(name)
    plt.xlabel('Number of threads')
    plt.ylabel('Speedup')

    # show plot
    plt.legend()
    plt.savefig("speedup.png")
    plt.show()


# main program
args = parse_cmd()
threads = read_binary(args.threads, str)
configs = read_binary(args.configs, str)
data = create_data(args.input, len(threads))
plot_speedup(data, threads, configs, args.diagram)
