#pragma once

#include "lib/id_queue.h"
#include "src/data-types.h"
#include "src/datastructure/timestamped_vector.hpp"
#include <vector>

/**
 * Implementation of Floyd-Warshall Algorithm.
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
 */
class FloydWarshall
{
  public:
    FloydWarshall(std::size_t size) : number_of_nodes(size), distance(size * size, INF), parent(size * size, INF)
    {
    }

    template <bool update_parents = false> void compute_all_distances(auto neighbors)
    {
        reset();
        NodeId n = number_of_nodes;

        // init distance entries
        for (NodeId i = 0; i < n; i++)
        {
            distance[get_index(i, i)] = 0;
            auto init_edge = [&](NodeId v, Distance weight) {
                distance[get_index(i, v)] = weight;
                parent[get_index(i, v)] = i;
            };
            neighbors(i, init_edge);
        }

        for (NodeId k = 0; k < n; k++)
        {
            for (NodeId i = 0; i < n; i++)
            {
                for (NodeId j = 0; j < n; j++)
                {
                    Distance alternative = distance[get_index(i, k)] + distance[get_index(k, j)];
                    if (alternative < distance[get_index(i, j)])
                    {
                        distance[get_index(i, j)] = alternative;
                        if constexpr (update_parents)
                            parent[get_index(i, j)] = parent[get_index(i, k)];
                    }
                }
            }
        }
    }

    // Returns distance between nodes.
    Distance get_distance(NodeId v, NodeId w)
    {
        return distance[get_index(v, w)];
    }

    // Returns parent of node, if v and w in same connected component, otherwise INF.
    NodeId get_parent(NodeId v, NodeId w)
    {
        return parent[get_index(v, w)];
    }

    // set number of nodes if you want to operate on a subgraph of smaller size. NodeIDs have to be mapped to 0, ..., n
    // - 1
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

    inline std::size_t get_index(int i, int j)
    {
        return i * number_of_nodes + j;
    }

    std::size_t number_of_nodes;
    TimestampedVector<Distance>
        distance; // represent 2d array as rolled out 1d array to use timestamped vector without modification
    TimestampedVector<NodeId> parent;
};
