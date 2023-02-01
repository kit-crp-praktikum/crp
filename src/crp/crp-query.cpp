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
        const int level = std::min(this->overlay->partition.find_level_differing(start, u),
                                   this->overlay->partition.find_level_differing(end, u));

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
                    auto dist = *overlay->get_distance(level, cellId, uId, vId);
                    relaxOp(neighbors[vId], dist);
                }
            }
        }
    };

    const auto &search_bwd = [&](NodeId u, auto relaxOp) {
        const int level = std::min(this->overlay->partition.find_level_differing(start, u),
                                   this->overlay->partition.find_level_differing(end, u));

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
                    relaxOp(neighbors[vId], *overlay->get_distanceT(level, cellId, vId, uId));
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
    for (unsigned i = 0; i < path.size() - 1; i++)
    {
        // consider u-v edge (possibly shortcut)
        NodeId u = path[i];
        NodeId v = path[i + 1];

        auto unpacked_path = _unpack(u, v);

        if (unpacked_path.size() > 2)
        {
            path.insert(path.begin() + i + 1, unpacked_path.begin() + 1, unpacked_path.end() - 1);
            // newly unpacked edge may be shortcut as well
            i--;
        }
    }

    return path;
}

std::vector<NodeId> CRPAlgorithm::_unpack(NodeId start, NodeId end)
{
    // fwd neighbors
    auto neighbors_in_cell_fwd = [&](NodeId v, auto f) {
        const int diff_start = partition.find_level_differing(start, v);
        const int diff_end = partition.find_level_differing(end, v);
        const int level = std::min(diff_start, diff_end);

        // Iterate graph edges
        for (auto [to, weight] : (*g)[v])
        {
            if (partition.find_level_differing(v, to) >= level) // in different cells on 'level'
            {
                f(to, weight);
            }
        }

        if (level >= 0)
        {
            // Iterate clique
            const CellId cellId_v = overlay->get_cell_for_node(v, level);
            const NodeId vId = overlay->get_internal_id(v, level);
            if (vId != g->num_nodes())
            {
                const std::span<NodeId> neighbors = overlay->get_border_nodes_for_cell(level, cellId_v);
                for (NodeId toId = 0; toId < neighbors.size(); toId++)
                {
                    if (toId != vId)
                    {
                        auto dist = *overlay->get_distance(level, cellId_v, vId, toId);
                        f(neighbors[toId], dist);
                    }
                }
            }
        }
    };
    // bwd neighbors
    auto neighbors_in_cell_bwd = [&](NodeId v, auto f) {
        const int level = std::min(partition.find_level_differing(start, v), partition.find_level_differing(end, v));

        // Iterate graph edges
        for (auto [to, weight] : (reverse)[v])
        {
            if (partition.find_level_differing(v, to) >= level)
            {
                f(to, weight);
            }
        }

        if (level >= 0)
        {
            // Iterate clique
            const CellId cellId_v = overlay->get_cell_for_node(v, level);
            const NodeId vId = overlay->get_internal_id(v, level);
            if (vId != g->num_nodes())
            {
                const std::span<NodeId> neighbors = overlay->get_border_nodes_for_cell(level, cellId_v);
                for (NodeId toId = 0; toId < neighbors.size(); toId++)
                {
                    if (toId != vId)
                    {
                        auto dist = *overlay->get_distance(level, cellId_v, toId, vId);
                        f(neighbors[toId], dist);
                    }
                }
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
