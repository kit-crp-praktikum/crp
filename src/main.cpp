#include "crp/crp.h"
#include "data-types.h"
#include "graph.h"
#include "lib/debug.h"
#include "lib/timer.h"
#include "lib/vector_io.h"
#include "partitioner/bfs-partitioner.h"
#include "partitioner/bipartitioner.h"
#include "partitioner/geo-data.h"
#include "partitioner/inertial_flow.hpp"
#include "partitioner/rec-partitioner.h"
#include <cstdlib>
#include <cstring>
#include <filesystem>

static void print_help()
{
    std::cout <<
        R"(Usage: crp [-i GRAPH_DIRECTORY] [-w EDGE_WEIGHTS] [-q QUERIES_DIR] [-l LEVELS] [-c CELLS_PER_LEVEL] [-p PARTITIONER]

Argument description:

-i, --input A directory containing the first_out, head, latitude, longitude and
                weights files
-w, --weight Can be either geo_distance or travel_time, depending on what we
                want to test.
-l, --levels The number of levels to use in CRP.
-c, --cells-per-level The number of cells in each recursive split.
-q, --queries A directory containing source and target files.
-v, --verify Instead of measuring the time for each query, verify that the
                answers are correct. The answers to the queries are assumed to
                be in the same folder as the queries, with name
                <weight_type>_length.
-p, --partitioner Which partitioner to use, currently bfs or inertial.
-h, --help Show a help message like this one.
)";
}

static int parse_integer_or_bail(char *param)
{
    char *end;
    int result = strtol(param, &end, 10);
    if (*end != '\0')
    {
        std::cout << "Invalid number given: \"" << param << "\"" << std::endl;
        std::exit(-1);
    }

    return result;
}

void check_file_exists_or_bail(std::string file)
{
    if (!std::filesystem::exists(file))
    {
        std::cout << "Missing graph data file " << file << std::endl;
        std::exit(-1);
    }
}

static void check_input_directory_valid(std::string data_dir, std::string weight, bool requires_geo_data)
{
    auto root = std::filesystem::path{data_dir};
    check_file_exists_or_bail(root / "first_out");
    check_file_exists_or_bail(root / "head");
    check_file_exists_or_bail(root / weight);
    if (requires_geo_data)
    {
        check_file_exists_or_bail(root / "latitude");
        check_file_exists_or_bail(root / "longitude");
    }
}

static void check_query_directory_exists(std::string data_dir)
{
    auto root = std::filesystem::path{data_dir};
    check_file_exists_or_bail(root / "source");
    check_file_exists_or_bail(root / "target");
}

static void check_verification_data_exists(std::string data_dir, std::string weight)
{
    auto root = std::filesystem::path{data_dir};
    check_file_exists_or_bail(root / (weight + "_length"));
}

int find_argument_index(int argc, char **argv, char short_arg, std::string long_arg)
{
    std::string sh = std::string("-") + short_arg;
    std::string ln = "--" + long_arg;

    for (int i = 1; i < argc; i++)
    {
        if (argv[i] == sh || argv[i] == ln)
        {
            return i;
        }
    }

    return -1;
}

// Find argument which MUST be provided
int find_required_argument(int argc, char **argv, char short_arg, std::string long_arg, bool requires_value)
{
    auto pos = find_argument_index(argc, argv, short_arg, long_arg);

    if (pos == -1)
    {
        std::cout << "Need to specify " << long_arg << std::endl;
        std::exit(-1);
    }

    if (requires_value && pos == argc - 1)
    {
        std::cout << "Missing value for argument " << long_arg << std::endl;
        std::exit(-1);
    }

    return pos;
}

auto inertial_flow_part = [](crp::Graph *g, partitioner::GeoData *geo_data, int nr_levels,
                             int nr_cells) -> crp::RecursivePartition {
    partitioner::InertialFlowPartitioner part(g->num_nodes(), 4, 0.25);
    partitioner::RecPartitioner rec(part, nr_cells, nr_levels);

    crp::RecursivePartition partition;
    partition.number_of_levels = nr_levels;
    partition.cells_per_level = nr_cells;

    auto list = g->to_list();
    partition.mask = rec.partition_rec(list, *geo_data);
    return partition;
};

auto bfs_part = [](crp::Graph *g, partitioner::GeoData *geo_data, int nr_levels,
                   int nr_cells) -> crp::RecursivePartition {
    partitioner::BfsPartitioner part;
    partitioner::RecPartitioner rec(part, nr_cells, nr_levels);

    crp::RecursivePartition partition;
    partition.number_of_levels = nr_levels;
    partition.cells_per_level = nr_cells;

    auto list = g->to_list();
    partition.mask = rec.partition_rec(list, *geo_data);
    return partition;
};

struct CmdLineParams
{
    crp::CRPAlgorithmParams algo_params;
    std::string data_dir;
    std::string query_dir;
    std::string weight_type;
    bool verify_query_results;
};

CmdLineParams load_parameters_from_cmdline(int argc, char **argv)
{
    // Check for help message
    if (find_argument_index(argc, argv, 'h', "help") != -1)
    {
        print_help();
        std::exit(0);
    }

    CmdLineParams params;

    // Check command-line arguments
    int pos;
    pos = find_required_argument(argc, argv, 'l', "levels", true);
    params.algo_params.number_of_levels = parse_integer_or_bail(argv[pos + 1]);

    pos = find_required_argument(argc, argv, 'c', "cells-per-level", true);
    params.algo_params.cells_per_level = parse_integer_or_bail(argv[pos + 1]);

    pos = find_required_argument(argc, argv, 'i', "input", true);
    params.data_dir = argv[pos + 1];
    pos = find_required_argument(argc, argv, 'w', "weight", true);
    params.weight_type = argv[pos + 1];
    check_input_directory_valid(params.data_dir, params.weight_type, true);

    pos = find_required_argument(argc, argv, 'q', "queries", true);
    params.query_dir = argv[pos + 1];
    check_query_directory_exists(params.query_dir);

    pos = find_argument_index(argc, argv, 'v', "verify");
    params.verify_query_results = (pos != -1);
    if (params.verify_query_results)
    {
        check_verification_data_exists(params.query_dir, params.weight_type);
    }

    pos = find_argument_index(argc, argv, 'p', "partitioner");

    std::string part = "inertial";
    if (pos != -1 && pos != argc - 1)
    {
        part = argv[pos + 1];
    }

    if (part == "intertial")
    {
        params.algo_params.partitioner = inertial_flow_part;
    }
    else
    {
        params.algo_params.partitioner = bfs_part;
    }

    params.algo_params.customizer = crp::customize_with_dijkstra;
    return params;
}

int main(int argc, char **argv)
{
    auto params = load_parameters_from_cmdline(argc, argv);

    // Load graph
    std::filesystem::path dir = params.data_dir;
    crp::Graph g(dir / "first_out", dir / "head", dir / params.weight_type);
    partitioner::GeoData geo_data(dir / "latitude", dir / "longitude");

    crp::CRPAlgorithm algorithm{params.algo_params};
    const uint64_t prepare_time = get_time_debug("preparation", [&] { algorithm.prepare(&g, &geo_data); });
    const uint64_t customization_duration = get_time_debug("customization", [&] { algorithm.customize(); });

    std::filesystem::path query_dir = params.query_dir;
    auto sources = load_vector<uint32_t>(query_dir / "source");
    auto targets = load_vector<uint32_t>(query_dir / "target");

    size_t nr_queries = sources.size();
    if (params.verify_query_results)
    {
        auto answers = load_vector<uint32_t>(query_dir / (params.weight_type + "_length"));
        size_t incorrect = 0;

        get_time_debug("queries", [&] {
            for (size_t i = 0; i < nr_queries; i++)
            {
                Distance answer = algorithm.query(sources[i], targets[i]);
                incorrect += answer != answers[i];
            }
        });

        std::cout << incorrect << " out of " << nr_queries << " are not correct." << std::endl;
    }
    else
    {
        std::cout << prepare_time << " " << customization_duration << std::endl;
        std::vector<uint64_t> query_times(nr_queries);

        get_time_debug("queries", [&] {
            for (size_t i = 0; i < nr_queries; i++)
            {
                query_times[i] = get_time([&] { algorithm.query(sources[i], targets[i]); });
            }
        });

        for (size_t i = 0; i < nr_queries; i++)
        {
            std::cout << query_times[i] << " ";
        }
    }
}
