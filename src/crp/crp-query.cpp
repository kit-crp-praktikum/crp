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
        const int level = std::min(partition.find_level_differing(start, u), partition.find_level_differing(end, u));

        // Iterate graph edges
        for (auto [v, dist] : (*g)[u])
        {
            if (partition.find_level_differing(u, v) >= level)
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
        const int level = std::min(partition.find_level_differing(start, u), partition.find_level_differing(end, u));

        // Iterate graph edges
        for (auto [v, dist] : reverse[u])
        {
            if (partition.find_level_differing(u, v) >= level)
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

    auto path_with_shortcuts = bidir_dijkstra->unpack(start, end, middle);
    std::vector<NodeId> path;
    path.push_back(start);

    for (unsigned i = 0; i < path_with_shortcuts.size() - 1; i++)
    {
        NodeId u = path_with_shortcuts[i];
        NodeId v = path_with_shortcuts[i + 1];
        // is u-v edge shortcut?
        bool is_shortcut = true;
        for (auto [to, w] : (*g)[u])
        {
            // edge u -> to with wieght w
            if (to == v)
            {
                is_shortcut = false;
                break;
            }
        }
        if (!is_shortcut)
        { // u-v real edge of path
            path.push_back(v);
        }
        else
        { // u-v shortcut, find shortest path with bidijkstra
            // search on lowest level where they are in the same cell, only inside the cell
            auto search_level = partition.find_level_differing(u, v) + 1;

            CellId cellId = overlay->get_cell_for_node(u, search_level);
            // fwd neighbors
            auto neighbors_in_cell_fwd = [&](NodeId v, auto f) {
                for (auto [to, weight] : (*g)[v])
                {
                    // all neighbors in g that have the same cellId
                    if (overlay->get_cell_for_node(to, search_level) == cellId)
                    {
                        f(to, weight);
                    }
                }
            };
            cellId = overlay->get_cell_for_node(v, search_level);
            // bwd neighbors
            auto neighbors_in_cell_bwd = [&](NodeId v, auto f) {
                for (auto [to, weight] : (reverse)[v])
                {
                    // all neighbors in g that have the same cellId
                    if (overlay->get_cell_for_node(to, search_level) == cellId)
                    {
                        f(to, weight);
                    }
                }
            };

            // run bidijkstra
            auto [u_v_middle, u_v_dist] =
                bidir_dijkstra->compute_distance_target<true>(u, v, neighbors_in_cell_fwd, neighbors_in_cell_bwd);
            auto u_v_path = bidir_dijkstra->unpack(u, v, u_v_middle);
            for (auto node : u_v_path)
            {
                if (node == u)
                    continue;
                path.push_back(node);
            }
        }
    }
    return path;
}
} // namespace crp
