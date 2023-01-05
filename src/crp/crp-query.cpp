#include "algorithms/dijkstra.hpp"
#include "crp/crp.h"
#include "data-types.h"

namespace crp
{
std::pair<NodeId, Distance> CRPAlgorithm::_query(NodeId start, NodeId end)
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

        // Iterate clique
        const CellId cellId = overlay->get_cell_for_node(u, level);
        const NodeId uId = overlay->get_internal_id(u, level);
        const std::span<NodeId> neighbors = overlay->get_border_nodes_for_cell(level, cellId);
        for (NodeId vId = 0; vId < neighbors.size(); vId++)
        {
            if (vId != uId)
            {
                relaxOp(neighbors[vId], *overlay->get_distance(level, cellId, uId, vId));
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

        // Iterate clique
        const CellId cellId = overlay->get_cell_for_node(u, level);
        const NodeId uId = overlay->get_internal_id(u, level);
        const std::span<NodeId> neighbors = overlay->get_border_nodes_for_cell(level, cellId);
        for (NodeId vId = 0; vId < neighbors.size(); vId++)
        {
            if (vId != uId)
            {
                relaxOp(neighbors[vId], *overlay->get_distance(level, cellId, vId, uId));
            }
        }
    };

    return bidir_dijkstra->compute_distance_target(start, end, search_fwd, search_bwd);
}

Distance CRPAlgorithm::query(NodeId start, NodeId end)
{
    return _query(start, end).second;
}

std::vector<NodeId> CRPAlgorithm::query_path(NodeId start, NodeId end, Distance &out_dist)
{
    auto [middle, distance] = _query(start, end);
    out_dist = distance;

    // TODO: implement path unpacking from bidir_dijkstra.
}
} // namespace crp
