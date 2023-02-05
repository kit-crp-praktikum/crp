#include "algorithms/dijkstra.hpp"
#include "crp/crp.h"
#include "data-types.h"
#include "graph.h"

namespace crp
{

auto CRPAlgorithm::get_fwd_scan(const NodeId &start, const NodeId &end)
{
    return [&](NodeId u, auto relaxOp) {
        const int level = std::min(this->overlay->partition.find_level_differing(start, u),
                                   this->overlay->partition.find_level_differing(end, u));

        // Iterate graph edges
        for (auto [v, dist] : fwd_remapped[u])
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
}

auto CRPAlgorithm::get_bwd_scan(const NodeId &start, const NodeId &end)
{
    return [&](NodeId u, auto relaxOp) {
        const int level = std::min(this->overlay->partition.find_level_differing(start, u),
                                   this->overlay->partition.find_level_differing(end, u));

        // Iterate graph edges
        for (auto [v, dist] : bwd_remapped[u])
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
}

// Query the shortest path from start to end using bidir_dijkstra.
// The result is the node where forward and backward search met and the computed shortest distance.
template <bool need_parents> std::pair<NodeId, Distance> CRPAlgorithm::_query(NodeId start, NodeId end)
{
    return bidir_dijkstra->compute_distance_target<need_parents>(start, end, get_fwd_scan(start, end),
                                                                 get_bwd_scan(start, end));
}

Distance CRPAlgorithm::query(NodeId start, NodeId end)
{
    start = node_mapping[start];
    end = node_mapping[end];
    return _query<false>(start, end).second;
}

std::vector<NodeId> CRPAlgorithm::query_path(NodeId start, NodeId end, Distance &out_dist)
{
    start = node_mapping[start];
    end = node_mapping[end];

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

    for (auto &x : path)
    {
        x = node_inverse_mapping[x];
    }

    return path;
}

std::vector<NodeId> CRPAlgorithm::_unpack(NodeId start, NodeId end)
{
    // run bidijkstra
    auto [middle, dist] =
        bidir_dijkstra->compute_distance_target<true>(start, end, get_fwd_scan(start, end), get_bwd_scan(start, end));
    auto unpacked_path = bidir_dijkstra->unpack(start, end, middle);
    return unpacked_path;
}
} // namespace crp
