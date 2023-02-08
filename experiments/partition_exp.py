import sys, os, argparse, math
from itertools import product

# command line interface
def parse_cmd():
    arg_parser = argparse.ArgumentParser()
    arg_parser.add_argument(      '--ex',                                                 required=True, help="path to executable of crp")
    arg_parser.add_argument(      '--run',                                                 required=True, help="file to put commands into")
    arg_parser.add_argument('-d', '--data',                                                 required=True, help="file to put dataframe into")
    arg_parser.add_argument('-o', '--out',                                                required=True  , help="directory to redirect std-out to")
    arg_parser.add_argument('-t', '--threads',         nargs="*",   type=int,  default=[1],  required=False, help="list of threads to use")
    arg_parser.add_argument(      '--metric',          nargs="*",             default=["travel_time"], required=False, help="list of metrics: travel_time, geo_distance")
    arg_parser.add_argument('-c', '--cells',           nargs="*", type=int,                 required=True,  help="list of cell")
    arg_parser.add_argument('-l', '--levels',          nargs="*", type=int,                 required=True,      help="list of level")
    arg_parser.add_argument('-m', '--min-cellsize',               type=int,  default=1000,   required=False,       help="remove configurations with num_nodes / cells^levels < min_cellsize")
    arg_parser.add_argument('-g', '--graphs',          nargs="*",                           required=True,       help="list paths to graph folder")
    arg_parser.add_argument('-p', '--partitioner',     nargs="*",                default=[],   required=False,       help="list partitioners to use: bfs, inertial, kahip")
    arg_parser.add_argument('-e', '--kahip-imbalance', nargs="*", type=float, default=[0.1],            required=False,       help="list imbalances for KaHIP")
    arg_parser.add_argument('--kahip-mode',      nargs="*",                   default=["strong"],           required=False,       help="list KaHIP modes: eco, ecosocial, fast, fastsocial, strong, strongsocial")
    arg_parser.add_argument('--iflow-size',      nargs="*", type=float,       default=[0.25],            required=False,       help="list of inertial flow ratios")
    arg_parser.add_argument('--iflow-lines',     nargs="*", type=int,         default=[4],              required=False,       help="list of inertial flow number of lines")
    arg_parser.add_argument('-M', '--max_level', type=bool,         default=False,              required=False,       help="for each cell size use maximal level")

    args = arg_parser.parse_args()
    print("args=%s" % args, file=sys.stderr)
    return args


config_param = ["graph", "metric", "cells", "levels", "threads", "partitioner", "kahip_imbalance", "kahip_mode", "iflow_size", "iflow_lines"]
temp_config = "./temp_config"
temp_time = "./temp_time"

def print_run(args, lines):
    with open(args.run, 'w') as f: 
        for line in lines:
            f.write(line + '\n')

def print_dataframe(list_configs, config_file, config_param):
    with open(config_file, 'w') as f: 
        for config in list_configs:
            line = [str(config[param]) if param in config else "NaN" for param in config_param]
            f.write(' '.join(line) + '\n')

def save_partition_time(args, lines, time_file):
    for i in range(len(lines)):
        lines[i] += " 2> temp && cat temp | grep partition_time | cut -d \"=\" -f2 >> " + time_file

def shared_param(ex, gr, w, t, lv, c):
    return f"{ex} -i {gr} -w {w} -t {t} -l {lv} -c {c}"

def get_graph_name(graph):
    return os.path.basename(os.path.normpath(graph))

def get_num_nodes(graph):
    # one node has 4 bytes
    bytes = os.stat(f"{graph}/first_out").st_size
    return bytes // 4

def graph_suffix(graph, metric, levels, cells):
    name = get_graph_name(graph)
    return f"{name}-{metric}-l{levels}-c{cells}"

def partition_name_bfs(graph, metric, levels, cells):
    return f"p-bfs-" + graph_suffix(graph, metric, levels, cells)

def partition_name_inertial(graph, metric, size, lines, levels, cells):
    return f"p-inertial-{size}-{lines}-" + graph_suffix(graph, metric, levels, cells)

def partition_name_kahip(graph, metric, mode, imbalance, levels, cells):
    return f"p-kahip-{mode}-{imbalance}-" +  graph_suffix(graph, metric, levels, cells)

def cells_too_small(args, graph, levels, cells):
    nodes = get_num_nodes(graph)
    cell_size = nodes // cells**levels
    return cell_size < args.min_cellsize

def shared_config(run, gr, w, t, lv, c):
    return {
            "run" : run,
            "graph": get_graph_name(gr),
            "metric": w,
            "threads" : t,
            "levels": lv,
            "cells": c
         }

# generates run file for bfs
def partition_bfs(args, lines, list_configs):
    for gr, w, t, lv, c in product(args.graphs, args.metric, args.threads, args.levels, args.cells):
        if cells_too_small(args, gr, lv, c):
            continue

        part_name = partition_name_bfs(gr, w, lv, c)
        line = shared_param(args.ex, gr, w, t, lv, c) + f" --dump-partition > {args.out}/{part_name}"
        lines.append(line)
        config = shared_config(1, gr, w, t, lv, c) 
        config["partitioner"] = "bfs"
        list_configs.append(config)

# generates run file for inertial flow
def partition_inertial(args, lines, list_configs):
    for gr, w, s, l, t, lv, c in product(args.graphs, args.metric, args.iflow_size, args.iflow_lines, args.threads, args.levels, args.cells):
        if cells_too_small(args, gr, lv, c):
            continue

        part_name = partition_name_inertial(gr, w, s, l, lv, c)
        line = shared_param(args.ex, gr, w, t, lv, c) + f" -p inertial --iflow-size {s} --iflow-lines {l} --dump-partition > {args.out}/{part_name}"
        lines.append(line)
        config = shared_config(1, gr, w, t, lv, c) 
        config["partitioner"] = "inertial"
        config["iflow_size"] = s
        config["iflow_lines"] = l
        list_configs.append(config)

# generates run file for kahip
def partition_kahip(args, lines, list_configs):
    for gr, w, mode, e, t, lv, c in product(args.graphs, args.metric, args.kahip_mode, args.kahip_imbalance, args.threads, args.levels, args.cells):
        if cells_too_small(args, gr, lv, c):
            continue
        part_name = partition_name_kahip(gr, w, mode, e, lv, c)
        line = shared_param(args.ex, gr, w, t, lv, c) + f" -p kahip --kahip-imbalance {e} --kahip-mode {mode} --dump-partition > {args.out}/{part_name}"
        lines.append(line)
        config = shared_config(1, gr, w, t, lv, c) 
        config["partitioner"] = "kahip"
        config["kahip_mode"] = mode
        config["kahip_imbalance"] = e
        list_configs.append(config)


def max_level(n, c):
	return math.floor(math.log(n) / math.log(c))

def max_level_partition_kahip(args, lines, list_configs):
    for gr, w, mode, e, t, c in product(args.graphs, args.metric, args.kahip_mode, args.kahip_imbalance, args.threads, args.cells):
        lv = max_level(get_num_nodes(gr) / args.min_cellsize, c)
        if lv <= 0:
            continue
        part_name = partition_name_kahip(gr, w, mode, e, lv, c)
        line = shared_param(args.ex, gr, w, t, lv, c) + f" -p kahip --kahip-imbalance {e} --kahip-mode {mode} --dump-partition > {args.out}/{part_name}"
        lines.append(line)
        config = shared_config(1, gr, w, t, lv, c) 
        config["partitioner"] = "kahip"
        config["kahip_mode"] = mode
        config["kahip_imbalance"] = e
        list_configs.append(config)



def partition_run(args):
    global temp_time, temp_config
    lines = list()
    list_configs = list()
    # generate commandline arguments for crp
    if args.max_level:
        if "kahip" in args.partitioner:
            max_level_partition_kahip(args, lines, list_configs)

    else:
        if "bfs" in args.partitioner:
            partition_bfs(args, lines, list_configs)
        if "inertial" in args.partitioner:
            partition_inertial(args, lines, list_configs)
        if "kahip" in args.partitioner:
            partition_kahip(args, lines, list_configs)

    save_partition_time(args, lines, temp_time)
    print_dataframe(list_configs, temp_config, config_param)

    # merge config and time file
    lines.append(f"paste {temp_config} {temp_time} > {args.data} ")
    # remove temporary files
    lines.append(f"rm {temp_config} {temp_time} temp")
    print_run(args, lines)

def main():
    args = parse_cmd()
    partition_run(args)

if __name__ == "__main__":
    main()
