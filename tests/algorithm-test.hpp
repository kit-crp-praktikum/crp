#pragma once

#include "crp/crp.h"
#include "data-types.h"
#include "graph-generator.hpp"
#include "graph.h"
#include "grid-graph.hpp"
#include "path-unpacker.h"
#include "shortest-path-algorithm.h"
#include <iostream>
#include <memory>
#define _ << " " <<
#define debug(x) #x << " = " << x

inline void setup_algorithm(std::unique_ptr<crp::CRPAlgorithmInterface> &algorithm, crp::Graph *g)
{
    algorithm->prepare(g);

    std::mt19937 mt(0);
    // Random integer in [l, r]
    const auto &rnd_integer = [&](int l, int r) { return std::uniform_int_distribution<int>(l, r)(mt); };
    for (auto &x : g->weights)
    {
        x = rnd_integer(1, 5);
    }

    algorithm->customize();
}

// NB: graph's weights are all 1, we generate real weights during customization
inline void test_algorithm(std::unique_ptr<crp::CRPAlgorithmInterface> algorithm,
                           crp::Graph g = generate_random_undirected_graph(100, 300, 1))
{
    setup_algorithm(algorithm, &g);

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

// NB: graph's weights are all 1, we generate real weights during customization
inline void test_algorithm_path(std::unique_ptr<crp::CRPAlgorithmInterface> algorithm,
                                crp::Graph g = generate_random_undirected_graph(100, 300, 1))
{
    setup_algorithm(algorithm, &g);

    for (int a = 0; a < g.num_nodes(); a++)
    {
        for (int b = 0; b < g.num_nodes(); b++)
        {
            Distance dist;
            auto path = algorithm->query_path(a, b, dist);

            if (dist == INF)
            {
                CHECK(path.empty());
            }
            else
            {
                CHECK(crp::isPathCorrect(&path, &g, dist) == crp::PathUnpackingResult::Ok);
            }
        }
    }
}
