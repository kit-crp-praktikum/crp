#pragma once

#include "data-types.h"
#include "node-slice.h"
#include <algorithm>
#include <optional>
#include <string>
#include <vector>
#include <iostream>

namespace crp
{
using AdjacencyList = std::vector<std::vector<std::pair<NodeId, Distance>>>;

class Graph
{
  public:
    // Empty graph
    Graph()
    {
    }

    // Load a graph from a set of files in the format of the praktikum.
    Graph(std::string first_out, std::string head, std::string weights);

    // Load a graph from an adjacency list, useful for simpler tests.
    Graph(const AdjacencyList &adj_list);

    detail::NodeSlice<Graph> operator[](size_t idx) const
    {
        return detail::NodeSlice(this, idx);
    }

    int num_nodes() const
    {
        return (int)first_out.size() - 1;
    }

    // Check whether the graph has a directed edge from u to v.
    // Uses linear search with the expectation that each node's degree is very small.
    inline std::optional<Distance> get_edge(NodeId u, NodeId v) const
    {
        auto begin = head.begin() + first_out[u];
        auto end = head.begin() + first_out[u + 1];
        auto it = std::find(begin, end, v);
        return it == end ? std::optional<Distance>{} : weights[it - head.begin()];
    }

    /**
     * Get a copy of the graph, but reversed.
     */
    Graph reversed() const;
  
    void clear() 
    {
        first_out.clear();
        head.clear();
        weights.clear();
    }

    // print graph for debugging purposes
    void print() 
    {
        for(NodeId v = 0; v < num_nodes(); v++)
        {
            std::cout << v << " -> ";
            for(NodeId e = first_out[v]; e < first_out[v + 1]; e++)
            {
                std::cout << "(" << head[e] << ", " << weights[e] << ") ";
            }
            std::cout << "\n";
        }
    }

  public:
    std::vector<uint32_t> first_out;
    std::vector<NodeId> head;
    std::vector<Distance> weights;
};
} // namespace crp
