#include "graph.h"
#include "partitioner/bfs-partitioner.h"
#include "partitioner/rec-partitioner.h"
#include <bitset>
#include <iostream>
#include <vector>

#include "data-types.h"
#include <iomanip>
#include <numeric>

// make private methods accessable for testing
#define private public
#include "crp/crp.h"
#include "tests/grid-graph.hpp"

int main()
{
    // setup begin -----------------------------------------------------------
    // copied graph and partition from os-test.cpp
    const int n = 8;
    crp::Graph g{generate_grid_graph(n)};

    // Generate simple 2-level partition
    crp::RecursivePartition rp;
    rp.cells_per_level = 4;
    rp.number_of_levels = 2;
    rp.mask.resize(g.num_nodes());

    dump_grid_graph_node_ids(n);
    rp.mask = generate_two_level_partition(n);

    auto partitioner = [&](crp::Graph *g, partitioner::GeoData *geodata, int number_of_levels, int cells_per_level) {
        return rp;
    };
    auto customizer = [&](crp::Graph *g, crp::OverlayStructure *) {
        for (uint32_t i = 0; i < g->weights.size(); i++)
        {
            g->weights[i] = 1;
        }
    };
    crp::CRPAlgorithmParams param = {rp.number_of_levels, rp.cells_per_level, partitioner, customizer};
    crp::CRPAlgorithm crp(param);
    partitioner::GeoData geo{};
    crp.prepare(&g, &geo);

    // print mutilevel partition
    std::cout << std::endl;
    for (int l = 0; l < param.number_of_levels; l++)
    {
        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < n; j++)
            {
                crp::CellId cell = crp.overlay->get_cell_for_node(encode(i, j), l);
                std::cout << std::setw(2) << std::setfill(' ') << cell << " ";
            }
            std::cout << std::endl;
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;

    crp.customize();

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

    Distance dist;
    std::cout << "run path query";
    auto path = crp.query_path(0, 46, dist);

    std::cout << "\n Unpacked path: ";
    for (auto node : path)
        std::cout << node << " ";
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
