#include "algorithms/dijkstra.hpp"
#include "crp/crp.h"
#include "crp/lru.hpp"
#include "data-types.h"
#include "graph.h"
#include <map>

namespace crp
{

inline int CRPAlgorithm::get_search_level(NodeId start, NodeId end, NodeId u)
{
    return std::min(this->overlay->partition.find_level_differing(start, u),
                    this->overlay->partition.find_level_differing(end, u));
}

auto CRPAlgorithm::get_fwd_scan(const NodeId &start, const NodeId &end)
{
    return [&](NodeId u, auto relaxOp) {
        const int level = get_search_level(start, end, u);

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
        const int level = get_search_level(start, end, u);

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

Path CRPAlgorithm::query_path(NodeId start, NodeId end, Distance &out_dist)
{
    return query_path_experimental(start, end, out_dist);
}

Path CRPAlgorithm::unpack_shortcut_one_level(NodeId u, NodeId v, LevelId level)
{
    const int big_cell = overlay->partition.find_cell_for_node(u, level);
    const auto &search_direction = [&](crp::Graph &directed_graph, auto &&get_distance) {
        return [&](NodeId x, auto relaxOp) {
            for (auto [y, w] : directed_graph[x])
            {
                const bool y_in_own_cell = overlay->partition.find_cell_for_node(y, level) == big_cell;
                const bool y_border_node =
                    (level <= 0 || (overlay->get_internal_id(y, level - 1) != directed_graph.num_nodes()));

                if (y_in_own_cell && y_border_node)
                {
                    relaxOp(y, w);
                }
            }

            if (level > 0)
            {
                const int lower_level = level - 1;
                const NodeId xId = overlay->get_internal_id(x, lower_level);
                if (xId != (NodeId)directed_graph.num_nodes())
                {
                    const CellId cell = overlay->get_cell_for_node(x, lower_level);
                    auto cell_nodes = overlay->get_border_nodes_for_cell(lower_level, cell);

                    for (NodeId yId = 0; yId < cell_nodes.size(); yId++)
                    {
                        relaxOp(cell_nodes[yId], get_distance(lower_level, cell, xId, yId));
                    }
                }
            }
        };
    };

    const auto &get_dist_fwd = [&](LevelId level, CellId cell, NodeId x, NodeId y) {
        return *overlay->get_distance(level, cell, x, y);
    };

    const auto &get_dist_bwd = [&](LevelId level, CellId cell, NodeId x, NodeId y) {
        return *overlay->get_distanceT(level, cell, y, x);
    };

    const auto &search_forward = search_direction(fwd_remapped, get_dist_fwd);
    const auto &search_backward = search_direction(bwd_remapped, get_dist_bwd);

    auto [middle, unused] = bidir_dijkstra->compute_distance_target<true>(u, v, search_forward, search_backward);
    auto path = bidir_dijkstra->unpack(u, v, middle);
    return path;
}

static LRUCache<std::pair<NodeId, NodeId>, Path> cache(1'000);

void CRPAlgorithm::set_cache_size(size_t size)
{
    if (size < 1)
    {
        std::cerr << "Cache size needs to be at least 1!" << std::endl;
        std::exit(-1);
    }

    cache = LRUCache<std::pair<NodeId, NodeId>, Path>{size};
}

void CRPAlgorithm::reset_cache_statistics()
{
    cache.reset_stats();
}

static bool is_cut_edge(crp::OverlayStructure &os, NodeId start, NodeId end)
{
    int unpack_level = os.partition.find_level_differing(start, end) + 1;
    if (unpack_level == os.partition.number_of_levels)
    {
        return true;
    }

    const NodeId invalid_iid = os.partition.mask.size();

    if (os.get_internal_id(start, unpack_level) == invalid_iid || os.get_internal_id(end, unpack_level) == invalid_iid)
    {
        return true;
    }

    return false;
}

// assumes it is not a cut edge!
static int determine_unpack_level(crp::OverlayStructure &os, NodeId start, NodeId end)
{
    int unpack_level = os.partition.find_level_differing(start, end) + 1;
    const NodeId invalid_iid = os.partition.mask.size();
    while (unpack_level + 1 < os.partition.number_of_levels &&
           os.get_internal_id(start, unpack_level + 1) != invalid_iid &&
           os.get_internal_id(end, unpack_level + 1) != invalid_iid)
    {
        ++unpack_level;
    }

    return unpack_level;
}

template <bool use_cache> void CRPAlgorithm::unpack_shortcut_recursive(NodeId u, NodeId v, LevelId level, Path &path)
{
    if (level < 0 || is_cut_edge(*overlay, u, v))
    {
        path.push_back(v);
        return;
    }

    if constexpr (use_cache)
    {
        auto cached = cache.get_value({u, v});
        if (cached)
        {
            path.insert(path.end(), cached->begin() + 1, cached->end());
            return;
        }
    }

    auto subpath = unpack_shortcut_one_level(u, v, level);
    Path unpacked;

    if (level == 0)
    {
        unpacked = subpath;
    }
    else
    {
        unpacked = {u};
        for (size_t i = 0; i + 1 < subpath.size(); i++)
        {
            unpack_shortcut_recursive<use_cache>(subpath[i], subpath[i + 1], level - 1, unpacked);
        }
    }

    path.insert(path.end(), unpacked.begin() + 1, unpacked.end());
    if constexpr (use_cache)
    {
        cache.push_value({u, v}, std::move(unpacked));
    }
}

template <bool use_cache> Path CRPAlgorithm::_query_path_original(NodeId start, NodeId end, Distance &out_dist)
{
    start = node_mapping[start];
    end = node_mapping[end];

    auto [meeting_point, dist] = _query<true>(start, end);
    out_dist = dist;

    auto [fwd_path, bwd_path] = bidir_dijkstra->unpack_separate(start, end, meeting_point);

    Path path = {start};
    for (size_t i = 0; i + 1 < fwd_path.size(); i++)
    {
        NodeId u = fwd_path[i];
        NodeId v = fwd_path[i + 1];
        LevelId level = get_search_level(start, end, u);
        unpack_shortcut_recursive<use_cache>(u, v, level, path);
    }

    for (size_t i = 0; i + 1 < bwd_path.size(); i++)
    {
        NodeId u = bwd_path[i];
        NodeId v = bwd_path[i + 1];
        LevelId level = get_search_level(start, end, v);
        unpack_shortcut_recursive<use_cache>(u, v, level, path);
    }

    for (auto &x : path)
    {
        x = node_inverse_mapping[x];
    }

    return path;
}

Path CRPAlgorithm::query_path_original(NodeId start, NodeId end, Distance &out_dist)
{
    return _query_path_original<false>(start, end, out_dist);
}

Path CRPAlgorithm::query_path_original_cache(NodeId start, NodeId end, Distance &out_dist)
{
    return _query_path_original<true>(start, end, out_dist);
}

std::vector<NodeId> CRPAlgorithm::query_path_experimental(NodeId start, NodeId end, Distance &out_dist)
{
    return _query_path_experimental<false>(start, end, out_dist);
}

std::vector<NodeId> CRPAlgorithm::query_path_experimental_cache(NodeId start, NodeId end, Distance &out_dist)
{
    return _query_path_experimental<true>(start, end, out_dist);
}

template <bool use_cache>
std::vector<NodeId> CRPAlgorithm::_query_path_experimental(NodeId start, NodeId end, Distance &out_dist)
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

        auto unpacked_path = _unpack<use_cache>(u, v);

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

template <bool use_cache> std::vector<NodeId> CRPAlgorithm::_unpack(NodeId start, NodeId end)
{
    if (is_cut_edge(*overlay, start, end))
    {
        return {start, end};
    }

    if constexpr (use_cache)
    {
        auto cached = cache.get_value({start, end});
        if (cached)
        {
            return *cached;
        }
    }

    // run bidijkstra
    const int unpack_level = determine_unpack_level(*overlay, start, end);
    const int unpack_cell = overlay->partition.find_cell_for_node(start, unpack_level);

    const auto &normal_fwd_scan = get_fwd_scan(start, end);
    const auto &normal_bwd_scan = get_bwd_scan(start, end);

    const auto &constrained_fwd_scan = [&](NodeId u, auto relaxOp) {
        normal_fwd_scan(u, [&](NodeId v, Distance dist) {
            if (overlay->get_cell_for_node(v, unpack_level) == unpack_cell)
            {
                relaxOp(v, dist);
            }
        });
    };

    const auto &constrained_bwd_scan = [&](NodeId u, auto relaxOp) {
        normal_bwd_scan(u, [&](NodeId v, Distance dist) {
            if (overlay->get_cell_for_node(v, unpack_level) == unpack_cell)
            {
                relaxOp(v, dist);
            }
        });
    };

    auto [middle, dist] =
        bidir_dijkstra->compute_distance_target<true>(start, end, constrained_fwd_scan, constrained_bwd_scan);

    auto unpacked_path = bidir_dijkstra->unpack(start, end, middle);

    if constexpr (use_cache)
    {
        auto cpy = unpacked_path;
        cache.push_value({start, end}, std::move(cpy));
    }

    return unpacked_path;
}
} // namespace crp
