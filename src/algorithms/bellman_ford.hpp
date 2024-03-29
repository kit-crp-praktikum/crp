#pragma once

#include "data-types.h"
#include "datastructure/timestamped_vector.hpp"
#include "lib/id_queue.h"
#include <iostream>
#include <vector>

/**
 * Implementation of Bellman-Ford Algorithm.
 * Takes as lambda a function that accepts a void(NodeId v, auto f) function.
 *
 * Example:
 * auto neighbors = [&](NodeId v, auto f) {
 *      for(auto [u, weight] : gr[v]) {
 *          f(u, weight);
 *      }
 *  };
 *  auto myf = [&](NodeId u, Distance weight) {
 *      std::cout << u << " " << weight << "\n";
 *  };
 *  neighbors(0, myf);
 *
 *
 * auto for_all_nodes = [&](auto f) {
 *     for(NodeId v = 0; v < n; v++) {
 *          f(v);
 *     }
 * };
 */
class BellmanFord
{
  public:
    BellmanFord(std::size_t size) : distance(size, INF), parent(size, INF), number_of_nodes(size)
    {
    }

    template <bool update_parents = false>
    void generic_compute_distance(NodeId start, auto neighbors, auto for_all_nodes)
    {
        reset();
        distance[start] = 0;
        bool changed = true;
        for (uint32_t i = 0; i < number_of_nodes - 1 && changed; i++)
        {
            changed = false;
            auto relax_neighbors = [&](NodeId v) {
                auto relax_operation = [&](NodeId u, Distance weight) {
                    Distance relaxed = distance[v] + weight;
                    if (weight < INF && relaxed < distance[u])
                    {
                        distance[u] = relaxed;
                        changed = true;
                        if constexpr (update_parents)
                            parent[u] = v;
                    }
                };
                neighbors(v, relax_operation);
            };
            for_all_nodes(relax_neighbors);
        }
    }

    template <bool update_parents = false> void compute_distance(NodeId start, auto neighbors)
    {
        auto for_all_nodes = [&](auto f) {
            for (NodeId v = 0; v < number_of_nodes; v++)
            {
                f(v);
            }
        };
        generic_compute_distance<update_parents>(start, neighbors, for_all_nodes);
    }

    // Return distance to node, if v is in same connected component as start node.
    // Otherwise INF.
    Distance tentative_distance(NodeId v)
    {
        return distance[v];
    }

    // Return parent of node, if v is in same connected component as start node.
    // Otherwise INF.
    NodeId get_parent(NodeId v)
    {
        return parent[v];
    }

    // set num nodes if you want to operate on a subgraph of smaller size. NodeIDs have to be mapped to 0, ..., n - 1
    void set_number_of_nodes(std::size_t n)
    {
        number_of_nodes = n;
    }

  private:
    void reset()
    {
        distance.reset();
        parent.reset();
    }

    TimestampedVector<Distance> distance;
    TimestampedVector<NodeId> parent;
    std::size_t number_of_nodes;
};
