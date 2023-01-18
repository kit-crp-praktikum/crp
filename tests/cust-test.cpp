
#include "data-types.h"
#include "graph.h"
#include "partitioner/geo-data.h"
#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <doctest/doctest.h>
#include <iomanip>
#include <numeric>

// make private methods accessable for testing
#define private public
#include "crp/crp.h"
#include "grid-graph.hpp"

// cannot pass customizer as function in doctest
void test_customization_on_grid_graph(int n, auto F, bool print = false)
{
    // n must be a power of two
    REQUIRE((n & (n - 1)) == 0);
    crp::Graph g{generate_grid_graph(n)};

    // Generate simple 2-level partition
    crp::RecursivePartition rp;
    rp.cells_per_level = 4;
    rp.number_of_levels = 2;
    rp.mask.resize(g.num_nodes());

    rp.mask = generate_two_level_partition(n);
    if (print)
        dump_grid_graph_node_ids(n);
    if (print)
        grid_graph_print_recursive_partition(n, rp);

    auto partitioner = [&](crp::Graph *g, partitioner::GeoData *, int number_of_levels, int cells_per_level) {
        return rp;
    };

    auto customizer = [&](crp::Graph *g, crp::OverlayStructure *overlay) { F(g, overlay); };

    crp::CRPAlgorithmParams param = {rp.number_of_levels, rp.cells_per_level, partitioner, customizer};
    crp::CRPAlgorithm crp(param);
    partitioner::GeoData geo_data{};
    crp.prepare(&g, &geo_data);
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
                    REQUIRE(d == manhatten_distance(v, u));
                }
            }
        }
    }
}

TEST_CASE("Grid graph n=8 customize_with_dijkstra")
{
    test_customization_on_grid_graph(8, crp::customize_with_dijkstra, true);
};

TEST_CASE("Grid graph n=8 customize_with_bellman_ford")
{
    test_customization_on_grid_graph(8, crp::customize_with_bellman_ford, 1);
};

TEST_CASE("Grid graph n=8 customize_dijkstra_rebuild")
{
    test_customization_on_grid_graph(8, crp::customize_dijkstra_rebuild, 2);
};

TEST_CASE("Grid graph n=8 customize_bellman_ford_rebuild")
{
    test_customization_on_grid_graph(8, crp::customize_bellman_ford_rebuild, 3);
};

TEST_CASE("Grid graph n=8 customize_floyd_warshall_rebuild")
{
    test_customization_on_grid_graph(8, crp::customize_floyd_warshall_rebuild, 4);
};
