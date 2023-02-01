#include "graph.h"
#include "lib/timer.h"
#include "lib/vector_io.h"
#include "partitioner/bfs-partitioner.h"
#include "partitioner/rec-partitioner.h"
#include "path-unpacker.h"
#include <bitset>
#include <iostream>
#include <vector>

#include "data-types.h"
#include <algorithm>
#include <fstream>
#include <iomanip>
#include <numeric>
#include <random>

// make private methods accessable for testing
#define private public
#include "crp/crp.h"
#include "tests/grid-graph.hpp"

int main_dijkstra_rank()
{
    auto path_first_out = "E:/Git/crp/data/karlsruhe/first_out";
    auto path_head = "E:/Git/crp/data/karlsruhe/head";
    auto path_weight = "E:/Git/crp/data/karlsruhe/travel_time";
    auto path_time_output = "E:/Git/crp/data/output/sample_times.txt";
    auto path_dist_output = "E:/Git/crp/data/output/sample_dists.txt";
    auto path_output = "E:/Git/crp/data/output/sample.txt";

    crp::Graph g(path_first_out, path_head, path_weight);
    int n = g.num_nodes();
    int num_start_nodes = 10;
    std::vector<std::vector<Distance>> dists(num_start_nodes);

    // Generate simple 2-level partition
    crp::RecursivePartition rp{2, 4};
    rp.mask.resize(g.num_nodes());
    std::vector<int> start_nodes(n);
    std::iota(start_nodes.begin(), start_nodes.end(), 0);   // fill the vector with [0,1,2,3...]
    random_shuffle(start_nodes.begin(), start_nodes.end()); // shuffle the vector
    start_nodes.resize(num_start_nodes);

    Dijkstra dij(n);
    auto neighbors = [&](NodeId v, auto f) {
        for (auto [u, weight] : g[v])
        {
            f(u, weight);
        }
    };

    int phantomlevels = 0;
    // crp::CRPAlgorithmParams param = {
    //     rp.number_of_levels, phantomlevels, rp.cells_per_level, partitioner, customizer};
    // crp::CRPAlgorithm crp(param);
    // partitioner::GeoData geo{};
    // crp.prepare(&g, &geo);

    // for (int i = 0; i < num_start_nodes; i++)
    //{
    //     std::cout<<"\nRound "<<i;
    //     dists[i].resize(n);
    //     int start = start_nodes[i];
    //     dij.compute_distance(start, neighbors);
    //     for (int node = 0; node < n; node++) {
    //         dists[i][node] = dij.tentative_distance(node);
    //     }
    // }

    // write distances to file
    std::fstream out(path_output, std::ios::out);
    std::fstream out_dist(path_dist_output, std::ios::out);
    std::fstream out_time(path_time_output, std::ios::out);

    out << n << " " << 1000 << std::endl;

    for (std::vector<Distance> dist_vec : dists)
    {
        for (Distance distance : dist_vec)
        {
            out_dist << distance << " ";
        }
        // file << std::endl;
    }

    // write  1000*n random times into file up to 10ms, normal distribution

    std::mt19937 rng(std::random_device{}());
    std::normal_distribution<float> normal_dis(0.4, 1.0); // mean 5 and std_dev 2

    for (int i = 0; i < num_start_nodes * n; i++)
    {
        float random_time = normal_dis(rng);
        if (random_time < 0 || random_time > 10)
        {
            i--;
            continue;
        } // discard values outside the range [0,10]
        out_time << random_time << " ";
    }

    out.close();
    out_dist.close();
    out_time.close();
}

int main()
{
    // setup begin --------------------------------------------------------------------
    // generate crp and graph
    const int n = 8;
    crp::CRPAlgorithmParams params;
    params.number_of_levels = 2;
    params.cells_per_level = 4;
    params.partitioner = [](crp::Graph *g, partitioner::GeoData *geo, int nr_levels, int cells_per_level) {
        crp::RecursivePartition rp{nr_levels, cells_per_level};
        rp.mask = generate_two_level_partition(n);
        return rp;
    };
    params.customizer = crp::customize_with_dijkstra;

    dump_grid_graph_node_ids(n);
    auto crp = std::make_unique<crp::CRPAlgorithm>(params);
    auto graph = crp::Graph{generate_grid_graph(8)};

    // test algo
    crp->prepare(&graph, {});

    std::mt19937 mt(0);
    // Random integer in [l, r]
    const auto &rnd_integer = [&](int l, int r) { return std::uniform_int_distribution<int>(l, r)(mt); };
    for (auto &x : graph.weights)
    {
        x = rnd_integer(1, 5);
        // x = 1;
    }
    //>>>>>>> a732b2a (path unpacking works on KA graph)
    //
    //    // print mutilevel partition
    //    std::cout << std::endl;
    //    for (int l = 0; l < params.number_of_levels; l++)
    //    {
    //        for (int i = 0; i < n; i++)
    //        {
    //            for (int j = 0; j < n; j++)
    //            {
    //                crp::CellId cell = crp->overlay->get_cell_for_node(encode(i, j), l);
    //                std::cout << std::setw(2) << std::setfill(' ') << cell << " ";
    //            }
    //            std::cout << std::endl;
    //        }
    //        std::cout << std::endl;
    //    }
    //    std::cout << std::endl;
    //
    //<<<<<<< HEAD
    //    crp.customize();
    //    // setup done ----------------------------------------------------------
    //||||||| parent of a732b2a (path unpacking works on KA graph)
    //    crp.customize();

    auto decode = [&](NodeId v) {
        std::pair<int, int> p = {v / n, v % n};
        return p;
    };

    auto manhatten_distance = [&](NodeId v, NodeId u) {
        auto [i, j] = decode(v);
        auto [x, y] = decode(u);
        return std::abs(i - x) + std::abs(j - y);
    };

    // setup done ----------------------------------------------------------
    //=======
    //
    //    crp->customize();
    //
    //    // finished setup --------------------------------------------------------
    //
    //    //print graph
    //    /*    for (int i = 0; i < graph.num_nodes(); i++) {
    //        std::cout<<"\n"<<i<<": ";
    //        for (auto [to, w] : graph[i]) { // edge i -> to with wieght w
    //            std::cout<<"{"<<to<<","<<w<<"} ";
    //        }
    //    }
    //    std::cout<<std::endl;
    //    */
    //>>>>>>> a732b2a (path unpacking works on KA graph)

    Distance dist;
    /*
        auto path = crp->query_path(36, 2, dist);

        std::cout<< 36 << " " << 2 << std::endl;
        std::cout<<"Dist: "<<dist<<std::endl;
        std::cout << "\n Unpacked path: ";
        for (auto node : path)
            std::cout << node << " ";
        std::cout<<std::endl;
        std::cout<<"@@@ path correctness: "<<(int)(crp::isPathCorrect(&path, &graph, dist))<<std::endl;
    */
    for (int a = 0; a < graph.num_nodes(); a++)
    {
        for (int b = 0; b < graph.num_nodes(); b++)
        {
            auto path = crp->query_path(a, b, dist);
            if (crp::isPathCorrect(&path, &graph, dist) != crp::PathUnpackingResult::Ok)
            {
                std::cout << a << " " << b << std::endl;
                std::cout << "Dist: " << dist << std::endl;
                std::cout << "\n Unpacked path: ";
                for (auto node : path)
                    std::cout << node << " ";
                std::cout << std::endl;
                std::cout << "@@@ path correctness: " << (int)(crp::isPathCorrect(&path, &graph, dist)) << std::endl;
            }
        }
    }
}

void main_partitioner()
{
    //    unsigned n = 32;
    //    std::vector<std::vector<std::pair<NodeId, Distance>>> gr(n);
    //    // path
    //    for (NodeId i = 0; i < n - 1; i++)
    //    {
    //        gr[i].push_back({i + 1, 1});
    //        gr[i + 1].push_back({i, 1});
    //    }
    //
    //    partitioner::BfsPartitioner bfs;
    //    partitioner::RecPartitioner recPart(&bfs);
    //
    //    std::vector<NodeId> *masks = recPart.partition_rec(3, 2, &gr);
    //
    //    for (unsigned i = 0; i < masks->size(); i++)
    //    {
    //        std::cout << "Node " << i << ": " << (std::bitset<9>)((*masks)[i]) << std::endl;
    //    }
}
