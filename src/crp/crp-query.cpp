#include "algorithms/dijkstra.hpp"
#include "crp/crp.h"
#include "data-types.h"
#include "graph.h"

namespace crp
{

// Query the shortest path from start to end using bidir_dijkstra.
// The result is the node where forward and backward search met and the computed shortest distance.
template <bool need_parents> std::pair<NodeId, Distance> CRPAlgorithm::_query(NodeId start, NodeId end)
{
    const auto &search_fwd = [&](NodeId u, auto relaxOp) {
        const int level = std::min(this->overlay->partition.find_level_differing(start, u), this->overlay->partition.find_level_differing(end, u));

        // Iterate graph edges
        for (auto [v, dist] : (*g)[u])
        {
            if (this->overlay->partition.find_level_differing(u, v) >= level)
            {
                relaxOp(v, dist);
            }
        }

        if (level >= 0)
        {
            // Iterate clique
            const CellId cellId = overlay->get_cell_for_node(u, level);
            const NodeId uId = overlay->get_internal_id(u, level);
            if (uId != g->num_nodes())
            {
                const std::span<NodeId> neighbors = overlay->get_border_nodes_for_cell(level, cellId);
                for (NodeId vId = 0; vId < neighbors.size(); vId++)
                {
                    if (vId != uId)
                    {
                        auto dist = *overlay->get_distance(level, cellId, uId, vId);
                        relaxOp(neighbors[vId], dist);
                    }
                }
            }
        }
    };

    const auto &search_bwd = [&](NodeId u, auto relaxOp) {
        const int level = std::min(this->overlay->partition.find_level_differing(start, u), this->overlay->partition.find_level_differing(end, u));

        // Iterate graph edges
        for (auto [v, dist] : reverse[u])
        {
            if (this->overlay->partition.find_level_differing(u, v) >= level)
            {
                relaxOp(v, dist);
            }
        }

        if (level >= 0)
        {
            // Iterate clique
            const CellId cellId = overlay->get_cell_for_node(u, level);
            const NodeId uId = overlay->get_internal_id(u, level);
            if (uId != g->num_nodes())
            {
                const std::span<NodeId> neighbors = overlay->get_border_nodes_for_cell(level, cellId);
                for (NodeId vId = 0; vId < neighbors.size(); vId++)
                {
                    if (vId != uId)
                    {
                        relaxOp(neighbors[vId], *overlay->get_distance(level, cellId, vId, uId));
                    }
                }
            }
        }
    };

    return bidir_dijkstra->compute_distance_target<need_parents>(start, end, search_fwd, search_bwd);
}

Distance CRPAlgorithm::query(NodeId start, NodeId end)
{
    return _query<false>(start, end).second;
}

std::vector<NodeId> CRPAlgorithm::query_path(NodeId start, NodeId end, Distance &out_dist)
{
    auto [middle, distance] = _query<true>(start, end);
    out_dist = distance;

    auto path = bidir_dijkstra->unpack(start, end, middle);

    // go through all levels, top to bottom
    // if two nodes are in the same cell on the current level, and in different ones one level below,
    // unpack the shortcut / edge
    for (int curr_level = this->overlay->partition.number_of_levels-1; curr_level >= 0; curr_level--) 
    {
        for (unsigned i = 0; i < path.size() - 1; i++)
        {
            // consider u-v edge (possibly shortcut)
            NodeId u = path[i];
            NodeId v = path[i + 1];

            int unpack_level = this->overlay->partition.find_level_differing(u,v);
            if (unpack_level == this->overlay->partition.number_of_levels-1) continue;  // edge is not shortcut
            if (unpack_level != curr_level-1) continue; 
            
            auto unpacked_path = _unpack(u, v, unpack_level);
            if (unpacked_path.size() > 2) 
            {
                // insert unpacked_path into path without u,v
                path.insert (path.begin() + i+1, unpacked_path.begin() +1, unpacked_path.end()-1);
                i += unpacked_path.size()-2;
            }
        }
    }
    return path;
}

std::vector<NodeId> CRPAlgorithm::_unpack (NodeId start, NodeId end, int level)
{
    // level is the utmost level where start,end in different cells
    // thus cellId is same for start, end on level+1
    CellId cellId = overlay->get_cell_for_node(start, level+1);
    // fwd neighbors
    auto neighbors_in_cell_fwd = [&](NodeId v, auto f) {
        for (auto [to, weight] : (*g)[v])
        {
            // all neighbors in g that are in the same cell one level above
            if (overlay->get_cell_for_node(to, level+1) == cellId)
            {
                f(to, weight);
            }
        }
    };
    // bwd neighbors
    auto neighbors_in_cell_bwd = [&](NodeId v, auto f) {
        for (auto [to, weight] : (reverse)[v])
        {
            // all neighbors in g that are in the same cell one level above
            if (overlay->get_cell_for_node(to, level+1) == cellId)
            {
                f(to, weight);
            }
        }
    };

    // run bidijkstra
    auto [middle, dist] =
        bidir_dijkstra->compute_distance_target<true>(start, end, neighbors_in_cell_fwd, neighbors_in_cell_bwd);
    auto unpacked_path = bidir_dijkstra->unpack(start, end, middle);
    return unpacked_path;
}
} // namespace crp
