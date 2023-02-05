import sys, os, argparse, copy
from itertools import product
from run_partition import shared_param, shared_config, print_dataframe, print_run
# command line interface
def parse_cmd():
    arg_parser = argparse.ArgumentParser()
    arg_parser.add_argument(      '--ex',                                                 required=True, help="path to executable of crp")
    arg_parser.add_argument(      '--run',                                                 required=True, help="file to put commands into")
    arg_parser.add_argument('-d', '--data',                                                 required=True, help="file to put dataframe into")
    arg_parser.add_argument('-o', '--out',                                                required=True  , help="directory to redirect std-out to")
    arg_parser.add_argument('-g',   '--graphs',                                             help="directory to graphs")
    arg_parser.add_argument('-p',  '--partitions',                                             help="directory to partitions used for customization")
    arg_parser.add_argument('-t', '--threads',         nargs="*",   type=int,  default=[1],  required=False, help="list of threads to use")
    arg_parser.add_argument('-r', '--runs',                  type=int, default=1,    required=False, help="number of repetitions")
    arg_parser.add_argument('-C', '--customizer',     nargs="*",                 default=[], required=False,       help="list customizers to use: dijkstra, bf, dijkstra-rebuild, bf-rebuild, fw-rebuild, bf-simd-rebuild")

    args = arg_parser.parse_args()
    print("args=%s" % args, file=sys.stderr)
    return args


config_param = ["run", "graph", "metric", "cells", "levels", "threads", "kahip_imbalance", "kahip_mode", "iflow_size", "iflow_lines", "customizer"]
temp_config = "./c_temp_config"
temp_time = "./c_temp_time"

def save_customizer_time(lines, time_file):
    for i in range(len(lines)):
        lines[i] += " 2> temp && cat temp | grep customizer_time | cut -d \"=\" -f2 >> " + time_file

# generates run file for customization
def customizer(args, lines, list_configs):
    partitions = [p for p in os.listdir(args.partitions) if os.path.isfile(os.path.join(args.partitions, p))]
    runs = [i for i in range(1, args.runs + 1)]
    for r, part, cust, t in product(runs, partitions, args.customizer, args.threads):
        parsed = part.split('-')
        (part_name, graph_name, w, lv, c) = (parsed[1], parsed[-4], parsed[-3], parsed[-2][1], parsed[-1][1])
        cust_name = list(copy.copy(parsed))
        cust_name[0] = 'c'
        cust_name = '-'.join(cust_name)
        gr = args.graphs + "/" + graph_name
        #dont need them actually -> todo fix dump-customization
        queries = args.graphs + "/" + graph_name + "/test"
        partition = args.partitions + "/" + part
        line = shared_param(args.ex, gr, w, t, lv, c) + f" -q {queries} -p {partition} -C {cust} --dump-customization > {args.out}/{cust_name}"
        lines.append(line)
        config = shared_config(1, gr, w, t, lv, c) 
        config["partitioner"] = part_name
        config["customizer"] = cust
        config["run"] = r
        if part_name == "kahip":
            config["kahip_mode"] = parsed[2]
            config["kahip_imbalance"] = parsed[3]
        elif part_name == "inertial":
            config["iflow_size"] = parsed[2]
            config["iflow_lines"] = parsed[3]
        list_configs.append(config)

def main():
    args = parse_cmd()
    lines = list()
    list_configs = list()

    customizer(args, lines, list_configs)
    save_customizer_time(lines, temp_time)
    print_dataframe(list_configs, temp_config, config_param)
    # merge config and time file
    lines.append(f"paste {temp_config} {temp_time} > {args.data} ")
    # remove temporary files
    lines.append(f"rm {temp_config} {temp_time} temp")
    print_run(args, lines)

if __name__ == "__main__":
    main()




