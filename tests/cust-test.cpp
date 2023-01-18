
#include "algorithms/dijkstra.hpp"
#include "data-types.h"
#include "graph.h"
#include "partitioner/geo-data.h"
#include "partitioner/inertial_flow.hpp"
#include "partitioner/rec-partitioner.h"
#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <doctest/doctest.h>
#include <iomanip>
#include <numeric>

// make private methods accessable for testing
#define private public
#include "crp/crp.h"
#include "graph-generator.hpp"
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

// cannot pass customizer as function in doctest
void test_customization_on_general_graph(int n, int m, auto F)
{
    auto g = generate_random_undirected_graph(n, m, 100);

    partitioner::GeoData geodata;

    geodata.longitude.resize(n);
    geodata.latitude.resize(n);

    std::mt19937 mt(0);
    for (int i = 0; i < n; i++)
    {
        geodata.longitude[i] = std::uniform_real_distribution<float>(0, 90)(mt);
        geodata.latitude[i] = std::uniform_real_distribution<float>(0, 90)(mt);
    }

    auto part = partitioner::InertialFlowPartitioner{(NodeId)n, 4, 0.25};
    auto rec_part = partitioner::RecPartitioner(part, 4, 4);

    crp::RecursivePartition rp;
    rp.cells_per_level = 4;
    rp.number_of_levels = 4;
    rp.mask = rec_part.partition_rec(g.to_list(), geodata);

    crp::OverlayStructure os(&g, rp);
    F(&g, &os);

    Dijkstra dijk(g.num_nodes());

    // check for each edge in overlay if distance is equal to manhatten distance in grid graph
    for (crp::LevelId level = 0; level < rp.number_of_levels; level++)
    {
        for (crp::CellId cell = 0; cell < rp.cells_per_level; cell++)
        {
            auto border_nodes = os.get_border_nodes_for_cell(level, cell);
            for (NodeId u = 0; u < border_nodes.size(); u++)
            {
                for (NodeId v = 0; v < border_nodes.size(); v++)
                {
                    auto d1 = *os.get_distance(level, cell, u, v);
                    dijk.compute_distance_target(border_nodes[u], border_nodes[v], [&](NodeId u, auto RelaxOp) {
                        for (auto [v, w] : g[u])
                        {
                            if (rp.find_cell_for_node(v, level) != cell)
                            {
                                continue;
                            }

                            RelaxOp(v, w);
                        }
                    });
                    auto d2 = dijk.tentative_distance(border_nodes[v]);
                    CHECK(d1 == d2);
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

TEST_CASE("Random graph with dijkstra customization")
{
    test_customization_on_general_graph(200, 600, crp::customize_with_dijkstra);
}

TEST_CASE("Random graph with dijkstra-rebuild customization")
{
    test_customization_on_general_graph(200, 600, crp::customize_dijkstra_rebuild);
}

TEST_CASE("Random graph with bf customization")
{
    test_customization_on_general_graph(200, 600, crp::customize_with_bellman_ford);
}

TEST_CASE("Random graph with bf-rebuild customization")
{
    test_customization_on_general_graph(200, 600, crp::customize_bellman_ford_rebuild);
}

TEST_CASE("Random graph with fw-rebuild customization")
{
    test_customization_on_general_graph(200, 600, crp::customize_floyd_warshall_rebuild);
}
