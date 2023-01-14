#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <doctest/doctest.h>

#include "algorithms/dijkstra.hpp"
#include "graph-generator.hpp"
#include "path-unpacker.h"

TEST_CASE("Path unpacking in Dijkstra")
{
    auto g = generate_random_undirected_graph(20, 50, 10);
    Dijkstra d(g.num_nodes());
    for (int i = 0; i < g.num_nodes(); i++)
    {
        for (int j = 0; j < g.num_nodes(); j++)
        {
            d.compute_distance_target<true>(i, j, [&](NodeId u, auto ForEachNeighbor) {
                for (auto [v, weight] : g[u])
                {
                    ForEachNeighbor(v, weight);
                }
            });

            auto path = d.unpack(i, j);
            auto dist = d.tentative_distance(j);
            CHECK(crp::isPathCorrect(&path, &g, dist) == crp::PathUnpackingResult::Ok);
        }
    }
}

TEST_CASE("Path unpacking in BiDirDijkstra")
{
    auto g = generate_random_undirected_graph(20, 50, 10);
    auto rev = g.reversed();

    BidirectionalDijstkra d(g.num_nodes());
    for (int i = 0; i < g.num_nodes(); i++)
    {
        for (int j = 0; j < g.num_nodes(); j++)
        {

            auto search_neighbors_fwd = [&](NodeId u, auto ForEachNeighbor) {
                for (auto [v, weight] : g[u])
                {
                    ForEachNeighbor(v, weight);
                }
            };

            auto search_neighbors_bwd = [&](NodeId u, auto ForEachNeighbor) {
                for (auto [v, weight] : g[u])
                {
                    ForEachNeighbor(v, weight);
                }
            };

            auto result = d.compute_distance_target<true>(i, j, search_neighbors_fwd, search_neighbors_bwd);

            auto path = d.unpack(i, j, result.first);
            CHECK(crp::isPathCorrect(&path, &g, result.second) == crp::PathUnpackingResult::Ok);
        }
    }
}
