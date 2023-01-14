#pragma once

#include "crp/crp.h"
#include "graph-generator.hpp"
#include "grid-graph.hpp"
#include <iostream>
#define _ << " " <<
#define debug(x) #x << " = " << x

inline void test_algorithm(std::unique_ptr<crp::CRPAlgorithmInterface> algorithm,
                           crp::Graph g = generate_random_undirected_graph(
                               100, 300,
                               1)) // NB: graph's weights are all 1, we generate real weights during customization
{
    algorithm->prepare(&g);

    std::mt19937 mt(0);
    // Random integer in [l, r]
    const auto &rnd_integer = [&](int l, int r) { return std::uniform_int_distribution<int>(l, r)(mt); };
    for (auto &x : g.weights)
    {
        x = rnd_integer(1, 5);
    }

    algorithm->customize();

    Dijkstra plain_dijkstra(g.num_nodes());
    const auto &shortest_path_dijkstra = [&](NodeId a, NodeId b) {
        plain_dijkstra.reset();
        plain_dijkstra.compute_distance_target(a, b, [&](NodeId u, auto ForEachNeighbor) {
            for (auto [v, weight] : g[u])
            {
                ForEachNeighbor(v, weight);
            }
        });

        return plain_dijkstra.tentative_distance(b);
    };

    for (int a = 0; a < g.num_nodes(); a++)
    {
        for (int b = 0; b < g.num_nodes(); b++)
        {
            REQUIRE(algorithm->query(a, b) == shortest_path_dijkstra(a, b));
        }
    }
}
