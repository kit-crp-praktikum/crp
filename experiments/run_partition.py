import sys, os, argparse
from itertools import product

# command line interface
def parse_cmd():
    arg_parser = argparse.ArgumentParser()
    arg_parser.add_argument(      '--ex',                                                 required=True, help="path to executable of crp")
    arg_parser.add_argument(      '--run',                                                 required=True, help="file to put commands into")
    arg_parser.add_argument(      '--config',                                                 required=True, help="file to put configs into")
    arg_parser.add_argument(      '--time',                                                 required=True, help="file to put times into")
    arg_parser.add_argument('-o', '--out',                                      default="."   , help="directory to print output to")
    arg_parser.add_argument('-t', '--threads',         nargs="*",   type=int,  default=[1],  required=False, help="list of threads to use")
    arg_parser.add_argument(      '--metric',          nargs="*",             default=["travel_time"], required=False, help="list of metrics: travel_time, geo_distance")
    arg_parser.add_argument('-r', '--runs',                  type=int, default=1,    required=False, help="number of repetitions")
    arg_parser.add_argument('-c', '--cells',           nargs="*", type=int,                 required=True,  help="list of cell")
    arg_parser.add_argument('-l', '--levels',          nargs="*", type=int,                 required=True,      help="list of level")
    arg_parser.add_argument('-m', '--min-cellsize',               type=int,  default=1000,   required=False,       help="remove configurations with num_nodes / cells^levels < min_cellsize")
    arg_parser.add_argument('-g', '--graphs',          nargs="*",                           required=True,       help="list paths to graph folder")
    arg_parser.add_argument('-p', '--partitioner',     nargs="*",                           required=True,       help="list partitioners to use: bfs, inertial, kahip")
    arg_parser.add_argument('-e', '--kahip-imbalance', nargs="*", type=float, default=[0.1],            required=False,       help="list imbalances for KaHIP")
    arg_parser.add_argument('--kahip-mode',      nargs="*",                   default=["strong"],           required=False,       help="list KaHIP modes: eco, ecosocial, fast, fastsocial, strong, strongsocial")
    arg_parser.add_argument('--iflow-size',      nargs="*", type=float,       default=[0.25],            required=False,       help="list of inertial flow ratios")
    arg_parser.add_argument('--iflow-lines',     nargs="*", type=int,         default=[4],              required=False,       help="list of inertial flow number of lines")

    args = arg_parser.parse_args()
    print("args=%s" % args, file=sys.stderr)
    return args


config_param = ["run", "graph", "metric", "cells", "levels", "threads", "partitioner", "kahip_imbalance", "kahip_mode", "iflow_size", "iflow_lines"]

def print_run(args, lines):
    with open(args.run, 'w') as f: 
        for line in lines:
            f.write(line + '\n')

def print_configs(args, list_configs):
    global config_param
    with open(args.config, 'w') as f: 
        for config in list_configs:
            line = [str(config[param]) if param in config else "NaN" for param in config_param]
            f.write(','.join(line) + '\n')

def save_partition_time(args, lines):
    for i in range(len(lines)):
        lines[i] += " 2> temp && cat temp | grep partition_time | cut -d \"=\" -f2 > " + args.time

def shared_param(ex, gr, w, t, lv, c):
    return f"{ex} -i {gr} -w {w} -t {t} -l {lv} -c {c}"

def get_graph_name(graph):
    return os.path.basename(os.path.normpath(graph))

def get_num_nodes(graph):
    # one node has 4 bytes
    bytes = os.stat(f"{graph}/first_out").st_size
    return bytes // 4

def partition_name_kahip(graph, metric, mode, imbalance, levels, cells):
    name = get_graph_name(graph)
    return f"p_kahip_{mode}_{imbalance}_{name}_{metric}_l{levels}_c{cells}"

def partition_name_inertial(graph, metric, size, lines, levels, cells):
    name = get_graph_name(graph)
    return f"p_inertial_{size}_{lines}_{name}_{metric}_l{levels}_c{cells}"

# generates run file for kahip
def partition_inertial(args, lines, list_configs):
    for gr, w, s, l, t, lv, c in product(args.graphs, args.metric, args.iflow_size, args.iflow_lines, args.threads, args.levels, args.cells):
        nodes = get_num_nodes(gr)
        cell_size = nodes // c**lv
        if cell_size < args.min_cellsize:
            continue

        part_name = partition_name_inertial(gr, w, s, l, lv, c)
        line = shared_param(args.ex, gr, w, t, lv, c) + f" -p inertial --iflow-size {s} --iflow-lines {l} --dump-partition > {args.out}/{part_name}"
        lines.append(line)
        config = {
            "run" : 1,
            "graph": get_graph_name(gr),
            "metric": w,
            "threads" : t,
            "partitioner":"inertial",
            "iflow_size":s,
            "iflow_lines": l,
            "levels": lv,
            "cells": c
         }
        list_configs.append(config)

def partition_kahip(args, lines, list_configs):
    for gr, w, mode, e, t, lv, c in product(args.graphs, args.metric, args.kahip_mode, args.kahip_imbalance, args.threads, args.levels, args.cells):
        nodes = get_num_nodes(gr)
        cell_size = nodes // c**lv
        if cell_size < args.min_cellsize:
            continue

        part_name = partition_name_kahip(gr, w, mode, e, lv, c)
        line = shared_param(args.ex, gr, w, t, lv, c) + f" -p kahip --kahip-imbalance {e} --kahip-mode {mode} --dump-partition > {args.out}/{part_name}"
        lines.append(line)
        config = {
            "run" : 1,
            "graph": get_graph_name(gr),
            "metric": w,
            "threads" : t,
            "partitioner":"kahip",
            "kahip_mode":mode,
            "kahip_imbalance": e,
            "levels": lv,
            "cells": c
         }
        list_configs.append(config)


def partition_run(args, lines, configs):
    if "kahip" in args.partitioner:
        partition_kahip(args, lines, configs)
    if "inertial" in args.partitioner:
        partition_inertial(args, lines, configs)

args = parse_cmd()
lines = list()
list_configs = list()


partition_run(args, lines, list_configs)
save_partition_time(args, lines)
print_run(args, lines)
print_configs(args, list_configs)

