#pragma once

#include <vector>
#include <string>
#include "data-types.h"
#include "node-slice.h"

namespace crp
{
class Graph
{
  public:
    // Load a graph from a set of files in the format of the praktikum.
    Graph(std::string first_out, std::string head, std::string weights);

    // Load a graph from an adjacency list, useful for simpler tests.
    Graph(const std::vector<std::vector<std::pair<NodeId, Distance>>>& adj_list);

    detail::NodeSlice<Graph> operator[] (size_t idx) const
    {
        return detail::NodeSlice(this, idx);
    }

    int num_nodes() const
    {
        return (int)first_out.size() - 1;
    }

  public:
    std::vector<uint32_t> first_out;
    std::vector<NodeId> head;
    std::vector<Distance> weights;
};
}
