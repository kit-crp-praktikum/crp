#pragma once

#include "geo-data.h"
#include "graph.h"
#include <cassert>

namespace partitioner
{
// Represents an induced subgraph of the full graph for CRP.
struct Subgraph
{
    crp::AdjacencyList graph;

    // A mapping from node IDs in the local graph to the node ID in the full graph.
    std::vector<NodeId> mapping;

    // The geodata for the graph, may be empty if the subgraph was created without GeoData
    GeoData geo_data;
};

// Create new (sub)subgraphs of the given graph.
//
// @param partition A vector containing an element for each node in @graph, indicating in which partition
//   it is.
// @param nr_parts The number of parts in partition.
template <class MaskType>
std::vector<Subgraph> generate_subgraphs(const Subgraph &graph, const std::vector<MaskType> &partition, int nr_parts)
{
    assert(graph.graph.size() == partition.size());

    const int n = graph.graph.size();
    std::vector<int> graph_sizes(nr_parts);
    std::vector<int> local_id(n);

    for (int i = 0; i < n; i++)
    {
        assert(partition[i] < nr_parts);
        local_id[i] = graph_sizes[partition[i]];
        graph_sizes[partition[i]]++;
    }

    std::vector<Subgraph> subgraphs(nr_parts);
    for (int i = 0; i < nr_parts; i++)
    {
        subgraphs[i].graph.resize(graph_sizes[i]);
        subgraphs[i].mapping.resize(graph_sizes[i]);

        if (!graph.geo_data.longitude.empty())
        {
            subgraphs[i].geo_data.latitude.resize(graph_sizes[i]);
            subgraphs[i].geo_data.longitude.resize(graph_sizes[i]);
        }
    }

    for (int x = 0; x < n; x++)
    {
        if (!graph.geo_data.longitude.empty())
        {
            subgraphs[partition[x]].geo_data.longitude[local_id[x]] = graph.geo_data.longitude[x];
            subgraphs[partition[x]].geo_data.latitude[local_id[x]] = graph.geo_data.latitude[x];
        }

        subgraphs[partition[x]].mapping[local_id[x]] = graph.mapping[x];

        for (auto [y, w] : graph.graph[x])
        {
            if (partition[x] == partition[y])
            {
                subgraphs[partition[x]].graph[local_id[x]].push_back({local_id[y], w});
            }
        }
    }

    return subgraphs;
}
} // namespace partitioner
