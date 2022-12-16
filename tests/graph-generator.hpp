#pragma once

#include "data-types.h"
#include "graph.h"
#include <algorithm>
#include <iterator>
#include <random>
#include <set>
#include <doctest/doctest.h>

/**
 * Generates a random connected undirected graph with N nodes and M edges.
 * Generated weights are between 1 and W.
 *
 * Note that N should be relatively small, because the method requires N^2 memory.
 */
inline crp::Graph generate_random_undirected_graph(NodeId N, int M, Distance W)
{
    REQUIRE(M <= N * (N - 1) / 2);

    std::random_device rd;
    std::mt19937 mt(rd());

    // Random integer in [l, r]
    const auto& rnd_integer = [&] (int l, int r) {
        return std::uniform_int_distribution<int>(l,r)(mt);
    };

    std::vector<std::vector<std::pair<NodeId, Distance>>> adjlist(N);
    const auto& add_edge = [&] (NodeId from, NodeId to, Distance weight) {
        adjlist[from].push_back({to, weight});
        adjlist[to].push_back({from, weight});
    };

    // A list of the edges of the full graph, we pick a subset of them.
    std::vector<std::pair<NodeId, NodeId>> available_edges;
    for (NodeId x = 0; x < N; x++) {
        for (NodeId y = x + 1; y < N; y++) {
            available_edges.push_back({x, y});
        }
    }


    // Step 1: to guarantee correctness, generate a tree with N vertices.
    std::set<std::pair<NodeId, NodeId>> tree_edges;
    for (NodeId x = 1; x < N; x++) {
        // Each node should get a parent

        NodeId p = rnd_integer(0, x - 1);
        Distance w = rnd_integer(1, W);
        add_edge(p, x, w);
        tree_edges.insert({p, x});
    }

    // Step 2: fill up the remaining edges
    std::shuffle(available_edges.begin(), available_edges.end(), mt);
    int count_edges = N - 1;
    while (count_edges < M) {
        auto [x, y] = available_edges.back();
        available_edges.pop_back();

        if (!tree_edges.count(available_edges.back())) {
            Distance w = rnd_integer(1, W);
            add_edge(x, y, w);
            ++count_edges;
        }
    }

    // Step 3: random shuffle the labels to mask the underlying tree structure
    std::vector<NodeId> labels(N);
    std::iota(labels.begin(), labels.end(), 0);
    std::shuffle(labels.begin(), labels.end(), mt);

    return crp::Graph{adjlist};
}
