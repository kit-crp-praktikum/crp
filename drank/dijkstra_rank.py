# This is a python script for generating a Dijkstra rank plot
# The 'times' input file contains the binary coded uint64 vector of query times
# Example: For ranksize=10 the first 10 queries are Rank 0 nodes, the next 10 are Rank 1, etc.


import matplotlib.pyplot as plt
import numpy as np
import argparse


# command line interface
def parse_cmd():
    arg_parser = argparse.ArgumentParser()
    arg_parser.add_argument('-r', '--ranksize', default=1, type=int, help='size of each rank')
    arg_parser.add_argument('-t', '--times', help="directory of file 'times' containing running times")
    arg_parser.add_argument('-n', '--name', default='Dijkstra Rank', help='diagram name')
    arg_parser.add_argument('-l', '--levels', type=int, default='0', help='number of levels')
    arg_parser.add_argument('-c', '--cellsperlevel', type=int, default='0', help='number of cells on each level')

    args = arg_parser.parse_args()
    print("args=%s" % args)
    return args


# helper method to read binary encoded files
def read_binary(path, type):
    with open(path, mode='rb') as file:
        numbers = np.fromfile(file, dtype=type)
    return numbers


# read times
def read_times(folder, type):
    path_times = folder + "/times"
    times = read_binary(path_times, type)
    return times


# create data for plot
def create_data(all_times, rank_size):
    data = []
    for i in range(0, len(all_times), rank_size):
        data.append(all_times[i:i + rank_size])
    return data


# make dijkstra rank box plot
# data contains lists of times for each rank
def plot_dijkstra_rank(data, name):
    rank_labels = ["2^" + str(i) for i in range(len(data))]

    # set plot size
    plt.figure(figsize=(16, 7))

    # convert y-axis to Logarithmic scale
    plt.yscale("log")
    plt.grid(which="minor", linestyle="--")

    # creating plot
    plt.boxplot(data, labels=rank_labels, patch_artist=True)

    plt.title(name)
    plt.xlabel('Dijkstra Rank')
    plt.ylabel('Query time (microsec)')

    # show plot
    plt.savefig("dijkstra-rank.png")
    plt.show()


# main program
args = parse_cmd()
all_times = read_times(args.times, np.uint64)
data = create_data(all_times, args.ranksize)
plot_dijkstra_rank(data, args.name)
