#pragma once

#include "graph.h"
#include "data-types.h"

namespace partitioner {
    using AdjacencyList = std::vector<std::vector<std::pair<NodeId, Distance>>>;
    /**
     * directed edges in both directions -> weight 2
     * directed edge in one direction -> weight 1
     * assumes node degree is small for most of the nodes
    */
    AdjacencyList make_undirected(crp::Graph &graph) {
        AdjacencyList new_graph(graph.num_nodes());
        auto add_edge = [&](NodeId v, NodeId w, Distance d) {
            new_graph[v].push_back({w, d});
            new_graph[w].push_back({v, d});
        };
        bool directed;
        for(NodeId v = 0; v < (uint32_t)graph.num_nodes(); v++) {
            for (auto [w, d1] : graph[v]) { //edge (v, w)
                directed = true;
                for (auto [u, d2] : graph[w]) {
                    if(v == u) {  //(u, w) == (v, w)?
                        directed = false;
                        break;
                    }
                } 
                add_edge(v, w, directed? 1 : 2);
            } 
        }
        return new_graph;
    }
}