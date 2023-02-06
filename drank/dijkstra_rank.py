# This is a python script for generating a Dijkstra rank plot
# The 'input' contains the list of binary coded uint64 vectors of query times, which are to be compared
# Example: For ranksize=10 the first 10 queries are Rank 0 nodes, the next 10 are Rank 1, etc.

# Example run configs:
# dijkstra_rank.py -r 10 -i {folder_path}/times_query {folder_path}/times_query_path -n query query_path


import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import pandas as pd
import argparse


# command line interface
def parse_cmd():
    arg_parser = argparse.ArgumentParser()
    arg_parser.add_argument('-r', '--ranksize', default=1, type=int, help='size of each rank')
    arg_parser.add_argument('-i', '--input', nargs="*", help="files containing running times to compare")
    arg_parser.add_argument('-n', '--name', nargs="*",  help="names of results being compared")
    arg_parser.add_argument('-d', '--diagram', default='Dijkstra Rank', help='diagram name')

    args = arg_parser.parse_args()
    print("args=%s" % args)
    return args


# helper method to read binary encoded files
def read_binary(path, type):
    with open(path, mode='rb') as file:
        numbers = np.fromfile(file, dtype=type)
    return numbers


# create data of one file
def create_single_dataframe(times, rank_size, name):
    data = []
    for i in range(len(times)):
        data.append([times[i], '2^' + str(i//rank_size), name])

    # convert list to DataFrame
    df = pd.DataFrame(data, columns=['time', 'rank', 'name'])
    return df


# create whole dataframe
def create_dataframe(input_files, rank_size, names):
    headers = ['time', 'rank', 'name']
    df = pd.DataFrame(columns=headers)
    for i in range(len(input_files)):
        times = read_binary(input_files[i], np.uint64)
        times_df = create_single_dataframe(times, rank_size, names[i])
        df = pd.concat([df, times_df], axis=0)
    return df


# make dijkstra rank box plot
# data contains lists of times for each rank
def plot_dijkstra_rank(data, name):
    # set plot size
    plt.figure(figsize=(16, 7))

    # convert y-axis to Logarithmic scale
    plt.yscale("log")
    plt.grid(which="minor", linestyle="--")

    # create a grouped boxplot
    sns.boxplot(x=data['rank'], y=data['time'], hue=data['name'])

    plt.title(name)
    plt.xlabel('Dijkstra Rank')
    plt.ylabel('Query time (microsec)')

    # show plot
    plt.savefig("dijkstra-rank.png")
    plt.show()


# main program
args = parse_cmd()
df = create_dataframe(args.input, args.ranksize, args.name)
plot_dijkstra_rank(df, args.diagram)
