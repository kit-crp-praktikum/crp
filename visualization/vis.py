import matplotlib.pyplot as plt
import numpy as np
import matplotlib.cm as cm
import geopandas
import pandas as pd
from shapely import wkt
from shapely.geometry import Point, LineString
import sys
import argparse
import colorsys

### IO ###
#command line interface
def parse_cmd():
	argParser = argparse.ArgumentParser()
	argParser.add_argument("-g", "--graph",               	        help="directory of graph and geodata")
	argParser.add_argument("-p", "--partition",           	        help="partition file")
	argParser.add_argument("-l", "--levels",   default=1, type=int, help="number of partition levels")
	argParser.add_argument("-c", "--cells",               type=int, help="number of cells per level")
	argParser.add_argument("-s", "--sample",   default=1, type=int, help="take only every nth datapoint for faster plotting")
	argParser.add_argument("-m", "--markersize",   default=4, type=int, help="size of the dots in the plot")
	argParser.add_argument("-x", "--max_level",   default=10**9, type=int, help="draws level 0,..., min(levels, max_levels) - 1")
	argParser.add_argument("-e", "--edge_list", help="binary file of edges to visualize, format: fwd_size bwd_size (edges of fwd) (edges of bwd)")
	argParser.add_argument("-w", "--linewidth", default=0.05, type=float, help="width of drawn edges")

	args = argParser.parse_args()
	print("args=%s" % args)
	return args

#helper method to read binary encoded files
def read_binary(path, type):
    with open(path, mode='rb') as file:
        numbers = np.fromfile(file, dtype=type)
    return numbers

#read longitude and latitude from graph folder
def read_coordinates(folder):
	path_longitude = folder + "/longitude"
	path_latitude = folder + "/latitude"
	longitude = read_binary(path_longitude, np.float32)
	latitude = read_binary(path_latitude, np.float32)
	return longitude, latitude

#read graph from first_out and head file
def read_graph(folder):
	path_first_out = folder + "/first_out"
	path_head = folder + "/head"
	first_out = read_binary(path_first_out, np.uint32)
	head = read_binary(path_head, np.uint32)

	#print n and m
	n = len(first_out) - 1
	m = len(head)
	print(f"graph with {n} nodes and {m} edges")

	#build adjaceny list
	graph = [list() for _ in range(n)]
	for v in range(n):
		for out in range(first_out[v], first_out[v + 1]):
			graph[v].append(head[out])
	return graph

def read_partition(path):
	return read_binary(path, np.uint32)
### IO ###

### Plotting ###
def get_colors(num_colors):
    colors=[]
    for i in np.arange(0., 360., 360. / num_colors):
        hue = i/360.
        lightness = (50 + np.random.rand() * 10)/100.
        saturation = (90 + np.random.rand() * 10)/100.
        colors.append(colorsys.hls_to_rgb(hue, lightness, saturation))
    return colors

#def plot_partition_and_border(longitude, latitude, partition, output_name, graph_name, markersize, graph):
def plot_partition_and_border(graph_data, plot_parameter, level):
	sample = plot_parameter.sample
	markersize = plot_parameter.markersize
	output_name = plot_parameter.output_name

	partition = graph_data.partition_at_level[level][::sample]
	longitude = graph_data.longitude[::sample]
	latitude = graph_data.latitude[::sample]

	#create colormap
	m = partition.max() + 1
	color_bag = get_colors(m)
	color = [color_bag[i] for i in partition]

	#dataframes for graph nodes
	df = pd.DataFrame({'Longitude':longitude, 'Latitude':latitude, 'Color':color })
	gdf = geopandas.GeoDataFrame(df, geometry=geopandas.points_from_xy(df.Longitude, df.Latitude))

	#dataframes for border nodes
	border_longitude, border_latiude = graph_data.border_node_coordinates(level)
	border_df = pd.DataFrame({'Longitude':border_longitude, 'Latitude':border_latiude, 'Color':'black' })
	border_gdf = geopandas.GeoDataFrame(border_df, geometry=geopandas.points_from_xy(border_df.Longitude, border_df.Latitude))

	#create map
	if graph_name == "germany":
		#draw borders of germany
		world = geopandas.read_file(geopandas.datasets.get_path('naturalearth_lowres'))
		ax = gdf.plot(color=gdf.Color, markersize=markersize, figsize=(24, 12))
		ax = border_gdf.plot(ax=ax, color=border_gdf.Color, markersize=markersize, figsize=(24, 12))
		world[world.name == 'Germany'].plot(ax=ax,  facecolor="none", edgecolor='black', linewidth=2)
	else:
		ax = gdf.plot(color=gdf.Color, markersize=markersize, figsize=(24, 12))
		ax = border_gdf.plot(ax=ax, color=border_gdf.Color, markersize=markersize, figsize=(24, 12))

	#save plot
	file_name = output_name + ".png"
	plt.savefig(file_name)
	print(f"--> wrote plot to {file_name}")
### Plotting ###

#converts multilevel partition into list of single level partitions 
def extract_multilevel_partition(partition, num_levels, cells_per_level):
	def cell_for_node(v, level):
	    bits_per_level = cells_per_level.bit_length() - 1
	    bits = (num_levels - level) * bits_per_level;
	    return partition[v] & ((1 << bits) - 1)
	partition_at_level = [np.asarray([cell_for_node(v, lv) for v in range(len(partition))]) for lv in range(num_levels)]
	return partition_at_level

class GraphData:
	def __init__(self, graph, name, partition_at_level, longitude, latitude): 
		self.graph = graph
		self.partition_at_level = partition_at_level
		self.name = name
		self.longitude = longitude 
		self.latitude = latitude
		
	def border_node_coordinates(self, level):
		partition = self.partition_at_level[level]
		border_x, border_y = list(), list()
		for v in range(len(self.graph)):
			for w in self.graph[v]:
				if partition[v] != partition[w]:
					border_x.append(self.longitude[v])
					border_y.append(self.latitude[v])
					break
		return border_x, border_y


class PlotParameter:
	def __init__(self, output_name, markersize, sample):
		self.output_name = output_name
		self.markersize = markersize
		self.sample = sample

# longitude x, latitude y
def draw_edges(edge_list, edge_colors, longitude, latitude, linewidth, output_name):
	assert(len(edge_list) == 2 * len(edge_colors))

	# write edge into LineString format
	edges= list()
	for i in range(0, len(edge_list), 2):
		v, w = edge_list[i], edge_list[i + 1]
		edges.append(LineString([(longitude[v], latitude[v]), (longitude[w], latitude[w])]))

	# create dataframe
	gdf = geopandas.GeoDataFrame({'Edges':edges,'Color':edge_colors }, geometry='Edges')

	# create plot
	world = geopandas.read_file(geopandas.datasets.get_path('naturalearth_lowres'))
	ax = gdf.plot(color=gdf.Color, figsize=(24, 12), linewidth=linewidth)
	world[world.name == 'Germany'].plot(ax=ax,  facecolor="none", edgecolor='black', linewidth=2)
	#linewidth=1
	file_name = output_name + ".png"
	plt.savefig(file_name)
	print(f"--> wrote plot to {file_name}")

def find_level(u, partition, number_of_levels, start, end):
	def find_level_differing(a, b):
	    diff = partition[a] ^ partition[b];
	    if not '1' in str(diff):
	    	return -1
	    else:
	    	first_diff = str(diff).index('1')
	    first_diff = first_diff // number_of_levels
	    return number_of_levels - first_diff - 1;
	return min(find_level_differing(start, u), find_level_differing(end, u))


args = parse_cmd()
graph = read_graph(args.graph)
graph_name = args.graph.split("/")[-1]
longitude, latitude = read_coordinates(args.graph)
multi_lv_partition = read_partition(args.partition)
partition_at_level = extract_multilevel_partition(multi_lv_partition, args.levels, args.cells)

assert(len(longitude) == len(latitude))
assert(len(longitude) == len(multi_lv_partition))

graph_data = GraphData(graph, graph_name, partition_at_level, longitude, latitude)

# visualize seach space
if args.edge_list != None:
	# (u, v) u <- v
	edge_list = read_binary(args.edge_list, np.uint32)
	num_fwd, num_bwd = edge_list[0], edge_list[1]
	edges = edge_list[2:]
	start = edges[1]
	end = edges[num_fwd + 1]

	edge_colors = ['red' for _ in range(num_fwd // 2)] + ['blue' for _ in range(num_bwd // 2)]
	output_name = "query_{}_lv_{}_c_{}".format(graph_data.name, args.levels, args.cells)
	draw_edges(edges, edge_colors, longitude, latitude, args.linewidth, output_name)

# visualize partition
else:
	levels = min(args.levels, args.max_level)
	for lv in range(levels):
		output_name = "{}_lv_{}_c_{}".format(graph_data.name, lv, args.cells)
		plot_parameters = PlotParameter(output_name, args.markersize, args.sample)
		plot_partition_and_border(graph_data, plot_parameters, lv)