
#include "data-types.h"
#include "graph.h"
#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <doctest/doctest.h>
#include <iomanip>
#include <numeric>

// make private methods accessable for testing
#define private public
#include "crp/crp.h"
#include "grid-graph.hpp"

TEST_CASE("Grid graph n=8")
{
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

    auto partitioner = [&](crp::Graph *g, int number_of_levels, int cells_per_level) { return rp; };
    auto customizer = crp::customize_with_dijkstra;

    crp::CRPAlgorithmParams param = {rp.number_of_levels, rp.cells_per_level, partitioner, customizer};
    crp::CRPAlgorithm crp(param);
    crp.prepare(&g);

    grid_graph_print_recursive_partition(n, rp);

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

    // check for each edge in overlay if distance is equal to manhatten distance in grid graph
    for (crp::LevelId level = 0; level < param.number_of_levels; level++)
    {
        for (crp::CellId cell = 0; cell < param.cells_per_level; cell++)
        {
            auto border_nodes = crp.overlay->get_border_nodes_for_cell(level, cell);
            for (NodeId v : border_nodes)
            {
                for (NodeId u : border_nodes)
                {
                    NodeId vId = crp.overlay->get_internal_id(v, level);
                    NodeId uId = crp.overlay->get_internal_id(u, level);
                    Distance d = *crp.overlay->get_distance(level, cell, vId, uId);
                    CHECK(d == manhatten_distance(v, u));
                }
            }
        }
    }
};
