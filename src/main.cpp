#include "cmdline-parsing.hpp"
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
#include "partitioner/kahip-wrapper.hpp"
#include "partitioner/preprocessing.hpp"
#include "partitioner/rec-partitioner.h"
#include "path-unpacker.h"
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <map>
#include <omp.h>
#include <random>

static void print_help()
{
    std::cout <<
        R"(Usage: crp [-i GRAPH_DIRECTORY] [-w EDGE_WEIGHTS] [-q QUERIES_DIR] [-l LEVELS] [-c CELLS_PER_LEVEL] [-p PARTITIONER] [-C CUSTOMIZER]

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
-o, --phantom_levels The number of phantom_levels to use in CRP.

-p, --partitioner Which partitioner to use, currently bfs or inertial. Default is inertial.
-C, --customizer Which customizer to use. Default is dijkstra, can be one of:
                     dijkstra, bf, dijkstra-rebuild, bf-rebuild, fw-rebuild, bf-simd-rebuild

-t, --threads The number of threads to use.
-r, --reorder-nodes reorder border nodes to achieve better customization and query time

-h, --help Show a help message like this one.
--dump-partition Only generate the partition of the graph and dump the data on stdout.
--dump-customization Dump customization data on stdout.
--path-unpacking <|original|original-cache|experimental|experimental-cache> Unpack paths during benchmark and
                verification instead of computing shortest distances only. If original is set, the experimental
                CRP algorithm is used.

--iflow-lines number of lines used for Inertial Flow
--iflow-size ratio in (0.0, 0.5) of block size of one part

--kahip-mode sets mode of KaHIPa: eco, ecosocial, fast, fastsocial, strong, strongsocial
--kahip-imbalance imbalance allowed during partitioning with KaHIPa

--warmup-queries N Number of random queries to execute to warm up cache.
--cache-size N The maximal number of paths to store in the path cache.
)";
}

static void check_query_directory_exists(std::string data_dir)
{
    auto root = std::filesystem::path{data_dir};
    check_file_exists_or_bail((root / "source").generic_string());
    check_file_exists_or_bail((root / "target").generic_string());
}

static void check_verification_data_exists(std::string data_dir, std::string weight)
{
    auto root = std::filesystem::path{data_dir};
    check_file_exists_or_bail((root / (weight + "_length")).generic_string());
}

auto inertial_flow_part = [](crp::Graph *g, partitioner::GeoData *geo_data, int nr_levels, int nr_cells,
                             partitioner::InertialFlowParameters params) -> crp::RecursivePartition {
    partitioner::InertialFlowPartitioner part(g->num_nodes(), params.number_of_lines, params.group_size);
    partitioner::RecPartitioner rec(part, nr_cells, nr_levels);

    crp::RecursivePartition partition{nr_levels, nr_cells};
    auto list = partitioner::make_undirected(g->to_list());
    partition.mask = rec.partition_rec(list, *geo_data);
    return partition;
};

auto bfs_part = [](crp::Graph *g, partitioner::GeoData *geo_data, int nr_levels,
                   int nr_cells) -> crp::RecursivePartition {
    partitioner::BfsPartitioner part;
    partitioner::RecPartitioner rec(part, nr_cells, nr_levels);

    crp::RecursivePartition partition{nr_levels, nr_cells};
    auto list = g->to_list();
    partition.mask = rec.partition_rec(list, *geo_data);
    return partition;
};

auto kahip_part = [](crp::Graph *g, partitioner::GeoData *geo_data, int nr_levels, int nr_cells,
                     partitioner::KaHIPParameters params) -> crp::RecursivePartition {
    crp::RecursivePartition partition{nr_levels, nr_cells};
    params.nr_levels = nr_levels;
    params.nr_cells_per_level = nr_cells;
    partition.mask = partitioner::kahip_partition_graph(g, params);
    return partition;
};

auto load_partition_from_file = [](std::string file, crp::Graph *g, partitioner::GeoData *geo_data, int nr_levels,
                                   int nr_cells) -> crp::RecursivePartition {
    crp::RecursivePartition partition{nr_levels, nr_cells};
    partition.mask = load_vector<uint32_t>(file);
    if (partition.mask.size() != (size_t)g->num_nodes())
    {
        std::cout << "Supplied partition has invalid size! (" << partition.mask.size() << " vs. " << g->num_nodes()
                  << " expected)" << std::endl;
        std::exit(-1);
    }

    return partition;
};

enum class OperationMode
{
    Benchmark,
    Verify,
    PartitionOnly,
    CustomizeOnly,
};

enum class PathUnpackingMode
{
    // Do not unpack paths at all.
    NoUnpacking,
    // Unpack paths with the experimental implementation from Nora.
    UnpackExperimental,
    // Unpack paths with the experimental implementation from Nora.
    UnpackExperimentalCache,
    // Unpack paths with the original algorithm from the CRP paper.
    UnpackOriginal,
    // Unpack paths with the original algorithm from the CRP paper with cache enabled.
    UnpackOriginalCache,
};

struct CmdLineParams
{
    crp::CRPAlgorithmParams algo_params;
    std::string data_dir;
    std::string query_dir;
    std::string weight_type;
    OperationMode mode = OperationMode::Benchmark;
    PathUnpackingMode unpack = PathUnpackingMode::NoUnpacking;
    int warmup_queries = 0;
    bool reorder_nodes = false;
};

static partitioner::KaHIPMode parse_kahip_mode(std::string value)
{
    static const std::map<std::string, partitioner::KaHIPMode> modes = {
        {
            "eco",
            partitioner::KaHIPMode::ECO,
        },
        {
            "ecosocial",
            partitioner::KaHIPMode::ECOSOCIAL,
        },
        {
            "fast",
            partitioner::KaHIPMode::FAST,
        },
        {
            "fastsocial",
            partitioner::KaHIPMode::FASTSOCIAL,
        },
        {
            "strong",
            partitioner::KaHIPMode::STRONG,
        },
        {
            "strongsocial",
            partitioner::KaHIPMode::STRONGSOCIAL,
        },
    };

    if (modes.count(value))
    {
        return modes.find(value)->second;
    }
    else
    {
        std::cerr << "Unrecognized kahip mode " << value << std::endl;
        std::exit(-1);
    }
}

void select_partitioner(int argc, char **argv, CmdLineParams &params)
{
    int pos = find_argument_index(argc, argv, 'p', "partitioner");
    std::string part = "inertial";
    if (pos != -1 && pos != argc - 1)
    {
        part = argv[pos + 1];
    }
    std::cerr << "partitioner=" << part << "\n";
    if (part == "inertial")
    {
        partitioner::InertialFlowParameters inertial;
        inertial.number_of_lines = 4;
        inertial.group_size = 0.25;

        int pos = find_argument_index(argc, argv, ' ', "iflow-lines");
        if (pos != -1 && pos != argc - 1)
        {
            inertial.number_of_lines = parse_integer_or_bail(argv[pos + 1]);
        }

        pos = find_argument_index(argc, argv, ' ', "iflow-size");
        if (pos != -1 && pos != argc - 1)
        {
            inertial.group_size = parse_double_or_bail(argv[pos + 1]);
        }
        std::cerr << "iflow-lines=" << inertial.number_of_lines << "\n";
        std::cerr << "iflow-group-size=" << inertial.group_size << "\n";
        using namespace std::placeholders;
        params.algo_params.partitioner = std::bind(inertial_flow_part, _1, _2, _3, _4, inertial);
    }
    else if (part == "bfs")
    {
        params.algo_params.partitioner = bfs_part;
    }
    else if (part == "kahip")
    {
        partitioner::KaHIPParameters kahip;

        int pos = find_argument_index(argc, argv, ' ', "kahip-imbalance");
        if (pos != -1 && pos != argc - 1)
        {
            kahip.imbalance = parse_double_or_bail(argv[pos + 1]);
        }
        else
        {
            kahip.imbalance = 0.1;
        }

        pos = find_argument_index(argc, argv, ' ', "kahip-mode");
        if (pos != -1 && pos != argc - 1)
        {
            kahip.mode = parse_kahip_mode(argv[pos + 1]);
            std::cerr << "kahip-mode=" << argv[pos + 1] << "\n";
        }
        else
        {
            kahip.mode = partitioner::KaHIPMode::STRONG;
            std::cerr << "kahip-mode=strong\n";
        }
        std::cerr << "kahip-imbalance=" << kahip.imbalance << "\n";
        using namespace std::placeholders;
        params.algo_params.partitioner = std::bind(kahip_part, _1, _2, _3, _4, kahip);
    }
    else
    {
        // Read partition from a file
        std::cerr << "Using partition from file " << part << std::endl;
        check_file_exists_or_bail(part);
        using namespace std::placeholders;
        params.algo_params.partitioner = std::bind(load_partition_from_file, part, _1, _2, _3, _4);
    }
}

static void read_customization_matrix_from_file(std::string filename, crp::Graph *g, crp::OverlayStructure *os)
{
    // Instead of reading the whole file at once, we read at most 1MB of data at once. This way, we avoid the
    // problem of having to store the data both in RAM.

    std::ifstream file_in(filename);

    std::vector<Distance> buffer(1e6);
    size_t pointer = 0;
    size_t cap = 0;

    const auto &next_distance_from_file = [&]() {
        if (pointer == cap)
        {
            file_in.read((char *)buffer.data(), sizeof(Distance) * buffer.size());
            cap = file_in.gcount() / sizeof(Distance);
            pointer = 0;
        }

        if (pointer == cap)
        {
            std::cout << "Customizer data file has incorrect size (insufficient data)!" << std::endl;
            std::exit(-1);
        }

        auto r = buffer[pointer];
        pointer++;
        return r;
    };

    for (crp::LevelId level = 0; level < os->get_number_of_levels(); level++)
    {
        for (crp::CellId cell = 0; cell < os->num_cells_in_level(level); cell++)
        {
            auto border_nodes = os->get_border_nodes_for_cell(level, cell);
            for (NodeId i = 0; i < border_nodes.size(); i++)
            {
                for (NodeId j = 0; j < border_nodes.size(); j++)
                {
                    *os->get_distance(level, cell, i, j) = next_distance_from_file();
                }
            }
        }
    }

    if (pointer != cap)
    {
        std::cout << "Customizer data file has incorrect size (extra data)!" << std::endl;
        std::exit(-1);
    }
}

static void dump_overlay_structure(crp::OverlayStructure *os)
{
    std::vector<Distance> data;
    for (crp::LevelId level = 0; level < os->get_number_of_levels(); level++)
    {
        for (crp::CellId cell = 0; cell < os->num_cells_in_level(level); cell++)
        {
            auto border_nodes = os->get_border_nodes_for_cell(level, cell);
            for (NodeId i = 0; i < border_nodes.size(); i++)
            {
                for (NodeId j = 0; j < border_nodes.size(); j++)
                {
                    data.push_back(*os->get_distance(level, cell, i, j));
                }
            }
        }
    }

    std::cout.write((char *)data.data(), sizeof(uint32_t) * data.size());
}

void select_customizer(int argc, char **argv, CmdLineParams &params)
{
    int pos = find_argument_index(argc, argv, 'C', "customizer");
    std::string customizer = "dijkstra-rebuild";
    if (pos != -1 && pos != argc - 1)
    {
        customizer = argv[pos + 1];
    }

    if (customizer == "dijkstra")
    {
        params.algo_params.customizer = crp::customize_with_dijkstra;
    }
    else if (customizer == "dijkstra-rebuild")
    {
        params.algo_params.customizer = crp::customize_dijkstra_rebuild;
    }
    else if (customizer == "bf")
    {
        params.algo_params.customizer = crp::customize_with_bellman_ford;
    }
    else if (customizer == "bf-rebuild")
    {
        params.algo_params.customizer = crp::customize_bellman_ford_rebuild;
    }
    else if (customizer == "fw-rebuild")
    {
        params.algo_params.customizer = crp::customize_floyd_warshall_rebuild;
    }
    else if (customizer == "bf-simd-rebuild")
    {
        params.algo_params.customizer = crp::customize_bf_simd_rebuild;
    }
    else
    {
        std::cerr << "Reading customizer data from file " << customizer << std::endl;
        check_file_exists_or_bail(customizer);
        using namespace std::placeholders;
        params.algo_params.customizer = std::bind(read_customization_matrix_from_file, customizer, _1, _2);
    }
    std::cerr << "customizer=" << customizer << "\n";
}

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

    pos = find_argument_index(argc, argv, 'o', "phantom_levels");
    params.algo_params.number_of_phantom_levels = 0;
    if (pos != -1)
    {
        params.algo_params.number_of_phantom_levels = parse_integer_or_bail(argv[pos + 1]);
    }

    pos = find_argument_index(argc, argv, 't', "threads");
    int number_of_threads = 1;
    if (pos != -1)
    {
        number_of_threads = parse_integer_or_bail(argv[pos + 1]);
    }

    // Explicitly disable dynamic teams
    omp_set_dynamic(0);
    omp_set_num_threads(number_of_threads);

    pos = find_required_argument(argc, argv, 'c', "cells-per-level", true);
    params.algo_params.cells_per_level = parse_integer_or_bail(argv[pos + 1]);

    pos = find_required_argument(argc, argv, 'i', "input", true);
    params.data_dir = argv[pos + 1];
    pos = find_required_argument(argc, argv, 'w', "weight", true);
    params.weight_type = argv[pos + 1];
    check_input_directory_valid(params.data_dir, params.weight_type, true);

    pos = find_argument_index(argc, argv, 'r', "reorder-nodes");
    if (pos != -1)
    {
        params.reorder_nodes = true;
    }

    std::string path = params.data_dir;
    if (path.back() == '/')
        path = path.substr(0, path.size() - 1);
    std::string graph_name = path.substr(path.find_last_of("/\\") + 1);
    std::cerr << "graph_name=" << graph_name << "\n";
    std::cerr << "weight=" << params.weight_type << "\n";
    std::cerr << "number_of_phantom_levels=" << params.algo_params.number_of_phantom_levels << "\n";
    std::cerr << "cells_per_level=" << params.algo_params.cells_per_level << "\n";
    std::cerr << "levels=" << params.algo_params.number_of_levels << "\n";
    std::cerr << "threads=" << number_of_threads << "\n";
    std::cerr << "reorder_nodes=" << params.reorder_nodes << "\n";

    select_partitioner(argc, argv, params);

    pos = find_argument_index(argc, argv, ' ', "dump-partition");
    if (pos != -1)
    {
        params.mode = OperationMode::PartitionOnly;
        return params;
    }

    select_customizer(argc, argv, params);
    pos = find_argument_index(argc, argv, ' ', "dump-customization");
    if (pos != -1)
    {
        params.mode = OperationMode::CustomizeOnly;
        return params;
    }

    pos = find_required_argument(argc, argv, 'q', "queries", true);
    params.query_dir = argv[pos + 1];
    check_query_directory_exists(params.query_dir);

    pos = find_argument_index(argc, argv, 'v', "verify");
    if (pos != -1)
    {
        params.mode = OperationMode::Verify;
        check_verification_data_exists(params.query_dir, params.weight_type);
    }

    pos = find_argument_index(argc, argv, ' ', "path-unpacking");
    if (pos != -1)
    {
        if (pos != argc - 1 && std::string("original") == argv[pos + 1])
        {
            std::cerr << "Using original path unpacking strategy" << std::endl;
            params.unpack = PathUnpackingMode::UnpackOriginal;
        }
        else if (pos != argc - 1 && std::string("original-cache") == argv[pos + 1])
        {
            std::cerr << "Using original path unpacking strategy with cache" << std::endl;
            params.unpack = PathUnpackingMode::UnpackOriginalCache;
        }
        else if (pos != argc - 1 && std::string("experimental-cache") == argv[pos + 1])
        {
            std::cerr << "Using experimental path unpacking strategy with cache" << std::endl;
            params.unpack = PathUnpackingMode::UnpackExperimentalCache;
        }
        else
        {
            std::cerr << "Using experimental path unpacking strategy" << std::endl;
            params.unpack = PathUnpackingMode::UnpackExperimental;
        }
    }

    pos = find_argument_index(argc, argv, ' ', "warmup");
    if (pos != -1 && pos != argc - 1)
    {
        params.warmup_queries = parse_integer_or_bail(argv[pos + 1]);
        std::cerr << "Using " << params.warmup_queries << " warmup queries." << std::endl;
    }

    pos = find_argument_index(argc, argv, ' ', "cache-size");
    if (pos != -1 && pos != argc - 1)
    {
        int cache_size = parse_integer_or_bail(argv[pos + 1]);
        std::cerr << "Setting cache size to " << cache_size << std::endl;
        crp::CRPAlgorithm::set_cache_size(cache_size);
    }

    return params;
}

static void handle_precompute_only(CmdLineParams &params, crp::Graph &g, partitioner::GeoData &geo_data)
{
    int total_number_of_levels = params.algo_params.number_of_levels + params.algo_params.number_of_phantom_levels;
    crp::RecursivePartition rp(0, 0);
    auto time_partitioner = get_time([&]() {
        rp = params.algo_params.partitioner(&g, &geo_data, total_number_of_levels, params.algo_params.cells_per_level);
    });
    {
        crp::OverlayStructure os(&g, rp);
        std::cerr << "\noverlay statistics\n";
        std::cerr << "total_border_nodes=" << os.total_border_nodes() << "\n";
        std::cerr << "total_memory_bytes=" << os.total_memory_bytes() << "\n";
        std::cerr << "largest_cell=" << os.largest_cell() << "\n\n";
    }
    std::cerr << "partition_time=" << time_partitioner << "\n";

    if (params.mode == OperationMode::CustomizeOnly)
    {
        std::vector<NodeId> map, revmap;
        if (params.reorder_nodes)
        {
            get_time_debug("node reordering", [&] { crp::CRPAlgorithm::reoder_nodes(g, rp, map, revmap); });
        }

        crp::OverlayStructure *os;

        auto time = get_time([&] {
            os = new crp::OverlayStructure(&g, rp);
            params.algo_params.customizer(&g, os);
            os->remove_phantom_levels(params.algo_params.number_of_phantom_levels);
        });

        std::cerr << "customizer_time=" << time << "\n";
        dump_overlay_structure(os);
    }
    else /* if (mode == OperationMode::PartitionOnly) */
    {
        // Dump partition data on stdout
        std::cout.write((char *)rp.mask.data(), sizeof(uint32_t) * rp.mask.size());
    }
}

static void run_warmup_queries(crp::CRPAlgorithm *algorithm, size_t n, CmdLineParams &params)
{
    if (params.warmup_queries == 0)
    {
        return;
    }

    std::random_device rd;
    std::mt19937 mt{rd()};
    auto dist = std::uniform_int_distribution<NodeId>(0, n - 1);

    get_time_debug("cache warmup queries", [&] {
        for (size_t i = 0; i < params.warmup_queries; i++)
        {
            Distance ans;
            size_t x = dist(mt);
            size_t y = dist(mt);

            switch (params.unpack)
            {
            case PathUnpackingMode::UnpackExperimentalCache:
                algorithm->query_path_experimental_cache(x, y, ans);
                break;
            case PathUnpackingMode::UnpackOriginalCache:
                algorithm->query_path_original_cache(x, y, ans);
                break;

            case PathUnpackingMode::UnpackExperimental:
            case PathUnpackingMode::NoUnpacking: // fallthrough
            case PathUnpackingMode::UnpackOriginal:
                std::cerr << "Path unpacking mode does not support warmup and cache!" << std::endl;
                std::exit(-1);
            }
        }
    });
}

static void run_queries_verify(crp::CRPAlgorithm *algorithm, crp::Graph *g, CmdLineParams &params)
{
    std::filesystem::path query_dir = params.query_dir;
    auto sources = load_vector<uint32_t>((query_dir / "source").generic_string());
    auto targets = load_vector<uint32_t>((query_dir / "target").generic_string());

    size_t nr_queries = sources.size();
    auto answers = load_vector<uint32_t>((query_dir / (params.weight_type + "_length")).generic_string());
    size_t correct = 0;

    nr_queries = 1000;

#ifdef PRINT_EDGES
    // we only want to have the path for one query
    nr_queries = 1;
#endif

    run_warmup_queries(algorithm, g->num_nodes(), params);
    const auto &run_query_and_check = [&](NodeId from, NodeId to, Distance expected) -> int // returns 1 on correct
    {
        Distance answer;
        Path path;
        // there is no valid path to check
        if(expected == INF && params.unpack != PathUnpackingMode::NoUnpacking) 
        {
            return 1;
        }

        switch (params.unpack)
        {
        case PathUnpackingMode::NoUnpacking:
            answer = algorithm->query(from, to);
            return answer == expected;
        case PathUnpackingMode::UnpackExperimental:
            path = algorithm->query_path_experimental(from, to, answer);
            break;
        case PathUnpackingMode::UnpackExperimentalCache:
            path = algorithm->query_path_experimental_cache(from, to, answer);
            break;
        case PathUnpackingMode::UnpackOriginal:
            path = algorithm->query_path_original(from, to, answer);
            break;
        case PathUnpackingMode::UnpackOriginalCache:
            path = algorithm->query_path_original_cache(from, to, answer);
            break;
        }

        auto check_result = crp::isPathCorrect(&path, g, answer);

        if (check_result == crp::PathUnpackingResult::EdgeMissing)
        {
            std::cout << "edge missing\n";
            return 0;
        }
        else if (check_result == crp::PathUnpackingResult::TotalLengthWrong)
        {
            std::cout << "total dist wrong\n";
            return 0;
        }
        else
        {
            return 1;
        }
    };

    get_time_debug("queries", [&] {
        for (size_t i = 0; i < nr_queries; i++)
        {
            correct += run_query_and_check(sources[i], targets[i], answers[i]);
        }
    });
#ifdef PRINT_EDGES
    return 0;
#endif
    std::cout << correct << " out of " << nr_queries << " queries are correct." << std::endl;
}

static void run_queries_benchmark(crp::CRPAlgorithm *algorithm, crp::Graph *g, CmdLineParams &params)
{
    std::filesystem::path query_dir = params.query_dir;
    auto sources = load_vector<uint32_t>((query_dir / "source").generic_string());
    auto targets = load_vector<uint32_t>((query_dir / "target").generic_string());
    size_t nr_queries = sources.size();

    std::vector<uint64_t> query_times(nr_queries);
    run_warmup_queries(algorithm, g->num_nodes(), params);
    get_time_debug("queries", [&] {
        Distance answer;
        for (size_t i = 0; i < nr_queries; i++)
        {
            query_times[i] = get_time([&] {
                switch (params.unpack)
                {
                case PathUnpackingMode::NoUnpacking:
                    algorithm->query(sources[i], targets[i]);
                    break;

                case PathUnpackingMode::UnpackExperimental:
                    algorithm->query_path_experimental(sources[i], targets[i], answer);
                    break;
                case PathUnpackingMode::UnpackExperimentalCache:
                    algorithm->query_path_experimental_cache(sources[i], targets[i], answer);
                    break;

                case PathUnpackingMode::UnpackOriginal:
                    algorithm->query_path_original(sources[i], targets[i], answer);
                    break;
                case PathUnpackingMode::UnpackOriginalCache:
                    algorithm->query_path_original_cache(sources[i], targets[i], answer);
                    break;
                }
            });
        }
    });

    std::cout.write((char *)query_times.data(), sizeof(uint64_t) * query_times.size());
}

int main(int argc, char **argv)
{
    std::cerr << "configuration\n";
    auto params = load_parameters_from_cmdline(argc, argv);
    std::cerr << "\n";
    // Load graph
    std::filesystem::path dir = params.data_dir;
    crp::Graph g((dir / "first_out").generic_string(), (dir / "head").generic_string(),
                 (dir / params.weight_type).generic_string());
    partitioner::GeoData geo_data((dir / "latitude").generic_string(), (dir / "longitude").generic_string());

    if (params.mode == OperationMode::PartitionOnly or params.mode == OperationMode::CustomizeOnly)
    {
        handle_precompute_only(params, g, geo_data);
        return 0;
    }

    crp::CRPAlgorithm algorithm{params.algo_params};
    get_time_debug("preparation", [&] { algorithm.prepare(&g, &geo_data); });
    get_time_debug("customization", [&] { algorithm.customize(params.reorder_nodes); });

    if (params.mode == OperationMode::Verify)
    {
        run_queries_verify(&algorithm, &g, params);
    }
    else
    {
        run_queries_benchmark(&algorithm, &g, params);
    }
}
