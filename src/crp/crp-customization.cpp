#include "algorithms/bellman_ford.hpp"
#include "algorithms/bellman_ford_simd.hpp"
#include "algorithms/dijkstra.hpp"
#include "algorithms/floyd_warshall.hpp"
#include "crp.h"
#include "graph.h"
#include "lib/timer.h"

#include <iostream>
#include <numeric>
#include <omp.h>

namespace crp
{

void CRPAlgorithm::reoder_nodes(crp::Graph &g, RecursivePartition &partition, std::vector<NodeId> &node_mapping,
                                std::vector<NodeId> &node_inverse_mapping)
{
    node_inverse_mapping.resize(g.num_nodes());
    std::iota(node_inverse_mapping.begin(), node_inverse_mapping.end(), 0);

    std::vector<int> border_for(g.num_nodes(), 0);
    for (NodeId i = 0; i < (NodeId)g.num_nodes(); i++)
    {
        for (auto [j, w] : g[i])
        {
            border_for[i] = std::max(border_for[i], partition.find_level_differing(i, j));
            border_for[j] = std::max(border_for[j], partition.find_level_differing(i, j));
        }
    }

    std::stable_sort(node_inverse_mapping.begin(), node_inverse_mapping.end(),
                     [&](NodeId a, NodeId b) { return border_for[a] > border_for[b]; });

    node_mapping.resize(g.num_nodes());
    for (NodeId i = 0; i < (NodeId)g.num_nodes(); i++)
    {
        node_mapping[node_inverse_mapping[i]] = i;
    }

    // Remap graph and partition

    crp::AdjacencyList fwd_list(g.num_nodes());
    RecursivePartitionMask remapped_mask(g.num_nodes());

    for (NodeId i = 0; i < (NodeId)g.num_nodes(); i++)
    {
        remapped_mask[node_mapping[i]] = partition.mask[i];
        for (auto [j, w] : g[i])
        {
            fwd_list[node_mapping[i]].push_back({node_mapping[j], w});
        }
    }

    g = crp::Graph{fwd_list};
    partition.mask = remapped_mask;
}

void CRPAlgorithm::customize()
{
    customize(false);
}

void CRPAlgorithm::customize(bool reorder_nodes)
{
    // Copy the original graph
    this->fwd_remapped = *g;
    if (reorder_nodes)
    {
        CRPAlgorithm::reoder_nodes(fwd_remapped, partition, node_mapping, node_inverse_mapping);
    }
    else
    {
        // Identity mapping
        node_mapping.resize(g->num_nodes());
        node_inverse_mapping.resize(g->num_nodes());
        std::iota(node_mapping.begin(), node_mapping.end(), 0);
        std::iota(node_inverse_mapping.begin(), node_inverse_mapping.end(), 0);
    }

    this->bwd_remapped = fwd_remapped.reversed();

    // This is the overlay on the recursive partition with phantom_levels.
    this->overlay = std::make_unique<OverlayStructure>(&this->fwd_remapped, this->partition);
    this->params.customizer(&this->fwd_remapped, overlay.get());
    this->overlay->remove_phantom_levels(this->params.number_of_phantom_levels);
    this->overlay->precompute_cliquesT();
}

uint32_t largest_cell_size(crp::OverlayStructure *overlay)
{
    uint32_t max_size = 0;
    // all nodes in level zero cell
    for (CellId cellId = 0; cellId < overlay->num_cells_in_level(0); cellId++)
    {
        max_size = std::max(max_size, (uint32_t)overlay->get_nodes_level0(cellId).size());
    }
    for (LevelId level = 1; level < overlay->get_number_of_levels(); level++)
    {
        for (CellId cellId = 0; cellId < overlay->num_cells_in_level(level); cellId++)
        {
            // union of border nodes of previous level
            uint32_t sum = 0;
            const std::span<CellId> child_cells = overlay->get_child_cellIds(level, cellId);
            for (CellId child_cell : child_cells)
            {
                sum += overlay->get_border_nodes_for_cell(level - 1, child_cell).size();
            }
            max_size = std::max(max_size, sum);
        }
    }
    return max_size;
}

void generic_customize(crp::Graph *g, crp::OverlayStructure *overlay, auto init_algo, auto compute_clique)
{
    int num_threads;
#pragma omp parallel
    {
        num_threads = omp_get_num_threads();
    }

    // init thread local datastructures
    uint32_t largest_cell = largest_cell_size(overlay);
    std::cerr << "-> largest cell in overlay: " << largest_cell << "\n";
    std::cerr << "customizer_time_per_level=";
    auto thread_algo = init_algo(num_threads, largest_cell);
    std::vector<crp::Graph> thread_graph(num_threads);
    std::vector<std::vector<NodeId>> thread_mapping(num_threads, std::vector<NodeId>(g->num_nodes()));

    // compute level 0 cliques from original graph
    LevelId level = 0;
    long long cur;
    cur = get_micro_time();
#pragma omp parallel for
    for (CellId cellId = 0; cellId < overlay->num_cells_in_level(level); cellId++)
    {
        int tid = omp_get_thread_num();
        auto neighbors_in_cell = [&](NodeId v, auto f) {
            for (auto [to, weight] : (*g)[v])
            {
                // all neighbors in g that have the same cellId
                if (overlay->get_cell_for_node(to, level) == cellId)
                {
                    f(to, weight);
                }
            }
        };
        compute_clique(level, cellId, neighbors_in_cell, thread_algo[tid], thread_graph[tid], thread_mapping[tid]);
    }
    cur = get_micro_time() - cur;
    std::cerr << cur;

    // compute level i cliques from level i - 1 cliques
    for (LevelId level = 1; level < overlay->get_number_of_levels(); level++)
    {
        cur = get_micro_time();
#pragma omp parallel for
        for (CellId cellId = 0; cellId < overlay->num_cells_in_level(level); cellId++)
        {
            int tid = omp_get_thread_num();
            auto neighbors_in_cell = [&](NodeId u, auto f) {
                const LevelId prev_level = level - 1;
                const CellId upper_cell_u = overlay->get_cell_for_node(u, level);
                const CellId lower_cell_u = overlay->get_cell_for_node(u, prev_level);

                // Iterate graph edges
                for (auto [to, weight] : (*g)[u])
                {
                    // nodes are border nodes in this level
                    if (upper_cell_u == overlay->get_cell_for_node(to, level) &&
                        lower_cell_u != overlay->get_cell_for_node(to, level - 1))
                    {
                        f(to, weight);
                    }
                }

                // Iterate clique
                const NodeId internal_u = overlay->get_internal_id(u, prev_level);
                const std::span<NodeId> neighbors = overlay->get_border_nodes_for_cell(prev_level, lower_cell_u);
                // invalid ID -> is no border node
                if (internal_u == (NodeId)g->num_nodes())
                    return;

                for (NodeId internal_to = 0; internal_to < neighbors.size(); internal_to++)
                {
                    if (internal_u != internal_to)
                    {
                        f(neighbors[internal_to],
                          *overlay->get_distance(prev_level, lower_cell_u, internal_u, internal_to));
                    }
                }
            };
            compute_clique(level, cellId, neighbors_in_cell, thread_algo[tid], thread_graph[tid], thread_mapping[tid]);
        }
        cur = get_micro_time() - cur;
        std::cerr << "," << cur;
    }
    std::cerr << "\n";
}

void customize_with_shortest_path(crp::Graph *g, crp::OverlayStructure *overlay, auto init_algo, auto one_to_all,
                                  auto distance)
{
    auto compute_clique = [&](LevelId level, CellId cellId, auto neighbors, auto sp_algo, crp::Graph,
                              std::vector<NodeId>) {
        const std::span<NodeId> border_nodes = overlay->get_border_nodes_for_cell(level, cellId);
        for (NodeId u : border_nodes)
        {
            one_to_all(sp_algo, u, neighbors, level, cellId);
            const NodeId internal_u = overlay->get_internal_id(u, level);
            for (NodeId to : border_nodes)
            {
                const NodeId internal_to = overlay->get_internal_id(to, level);
                *overlay->get_distance(level, cellId, internal_u, internal_to) = distance(sp_algo, to);
            }
        }
    };
    generic_customize(g, overlay, init_algo, compute_clique);
}

void customize_with_dijkstra(crp::Graph *g, crp::OverlayStructure *overlay)
{
    auto init_algo = [&](int num_threads, uint32_t) {
        return std::vector<Dijkstra>(num_threads, Dijkstra(g->num_nodes()));
    };
    auto one_to_all = [&](Dijkstra &algo, NodeId u, auto neighbors, LevelId, CellId) {
        algo.compute_distance<false>(u, neighbors);
    };
    auto distance = [&](Dijkstra &algo, NodeId to) { return algo.tentative_distance(to); };

    customize_with_shortest_path(g, overlay, init_algo, one_to_all, distance);
}

void customize_with_bellman_ford(crp::Graph *g, crp::OverlayStructure *overlay)
{
    auto for_all_nodes = [&](LevelId level, CellId cellId, auto f) {
        if (level == 0)
        {
            // all nodes in level zero cell
            const std::span<NodeId> nodes_in_cell = overlay->get_nodes_level0(cellId);
            for (NodeId v_orig : nodes_in_cell)
            {
                f(v_orig);
            }
        }
        else
        {
            // union of border nodes of previous level
            const std::span<CellId> child_cells = overlay->get_child_cellIds(level, cellId);
            for (CellId child_cell : child_cells)
            {
                const std::span<NodeId> border_child = overlay->get_border_nodes_for_cell(level - 1, child_cell);
                for (NodeId v_orig : border_child)
                {
                    f(v_orig);
                }
            }
        }
    };
    // we don't set bellman ford size to subgraph size, but it will probably converge earlier -> worst case 1 iteration
    // more than necessary
    auto init_algo = [&](int num_threads, uint32_t) {
        return std::vector<BellmanFord>(num_threads, BellmanFord(g->num_nodes()));
    };
    auto one_to_all = [&](BellmanFord &algo, NodeId u, auto neighbors, LevelId level, CellId cellId) {
        auto nodes = [&](auto f) { for_all_nodes(level, cellId, f); };
        algo.generic_compute_distance<false>(u, neighbors, nodes);
    };
    auto distance = [&](BellmanFord &algo, NodeId to) { return algo.tentative_distance(to); };
    customize_with_shortest_path(g, overlay, init_algo, one_to_all, distance);
}

void customize_rebuild_graph(crp::Graph *g, crp::OverlayStructure *overlay, auto init_algo, auto compute_clique_entries)
{
    auto for_all_nodes = [&](LevelId level, CellId cellId, auto f) {
        if (level == 0)
        {
            // all nodes in level zero cell
            const std::span<NodeId> nodes_in_cell = overlay->get_nodes_level0(cellId);
            for (NodeId v_orig : nodes_in_cell)
            {
                f(v_orig);
            }
        }
        else
        {
            // union of border nodes of previous level
            const std::span<CellId> child_cells = overlay->get_child_cellIds(level, cellId);
            for (CellId child_cell : child_cells)
            {
                const std::span<NodeId> border_child = overlay->get_border_nodes_for_cell(level - 1, child_cell);
                for (NodeId v_orig : border_child)
                {
                    f(v_orig);
                }
            }
        }
    };

    auto build_graph = [&](LevelId level, CellId cellId, auto neighbors_in_cell, crp::Graph &cell_graph,
                           std::vector<NodeId> &original_to_local_node) {
        cell_graph.clear();

        // compute mapping between cell_graph and original graph
        NodeId local_id = 0;
        auto compute_mapping = [&](NodeId v_orig) {
            original_to_local_node[v_orig] = local_id;
            local_id++;
        };
        for_all_nodes(level, cellId, compute_mapping);

        // build up adjacency-array
        auto add_neighbors = [&](NodeId v_orig) {
            cell_graph.first_out.push_back(cell_graph.head.size());
            auto add_edge = [&](NodeId to, Distance weight) {
                cell_graph.head.push_back(original_to_local_node[to]);
                cell_graph.weights.push_back(weight);
            };
            neighbors_in_cell(v_orig, add_edge);
        };
        for_all_nodes(level, cellId, add_neighbors);
        // Sentinel at the end
        cell_graph.first_out.push_back(cell_graph.head.size());
    };

    auto compute_clique = [&](LevelId level, CellId cellId, auto neighbors_in_cell, auto sp_algo,
                              crp::Graph &cell_graph, std::vector<NodeId> &original_to_local_node) {
        build_graph(level, cellId, neighbors_in_cell, cell_graph, original_to_local_node);
        auto neighbor_in_cell_graph = [&](NodeId v, auto f) {
            for (auto [u, weight] : cell_graph[v])
            {
                f(u, weight);
            }
        };
        auto map_to_local = [&](NodeId v) { return original_to_local_node[v]; };
        compute_clique_entries(level, cellId, cell_graph.num_nodes(), neighbor_in_cell_graph, map_to_local, sp_algo);
    };
    generic_customize(g, overlay, init_algo, compute_clique);
}

void customize_shortest_path_rebuild(crp::Graph *g, crp::OverlayStructure *overlay, auto init_algo, auto init_size,
                                     auto one_to_all, auto distance)
{
    auto compute_clique_entries = [&](LevelId level, CellId cellId, NodeId subgraph_size, auto neighbors,
                                      auto map_to_local, auto sp_algo) {
        init_size(sp_algo, subgraph_size, neighbors);
        const std::span<NodeId> border_nodes = overlay->get_border_nodes_for_cell(level, cellId);
        for (NodeId u : border_nodes)
        {
            const NodeId internal_u = map_to_local(u);
            one_to_all(sp_algo, internal_u, neighbors);
            for (NodeId to : border_nodes)
            {
                const NodeId internal_to = map_to_local(to);
                const NodeId clique_u = overlay->get_internal_id(u, level);
                const NodeId clique_to = overlay->get_internal_id(to, level);
                *overlay->get_distance(level, cellId, clique_u, clique_to) = distance(sp_algo, internal_u, internal_to);
            }
        }
    };
    customize_rebuild_graph(g, overlay, init_algo, compute_clique_entries);
}

void customize_dijkstra_rebuild(crp::Graph *g, crp::OverlayStructure *overlay)
{
    auto init_algo = [&](int num_threads, uint32_t largest_cell) {
        return std::vector<Dijkstra>(num_threads, Dijkstra(largest_cell));
    };
    auto init_size = [&](Dijkstra, NodeId, auto) { return; };
    auto one_to_all = [&](Dijkstra &algo, NodeId u, auto neighbors) { algo.compute_distance<false>(u, neighbors); };
    auto distance = [&](Dijkstra &algo, NodeId, NodeId to) { return algo.tentative_distance(to); };

    customize_shortest_path_rebuild(g, overlay, init_algo, init_size, one_to_all, distance);
}

// same code as dijkstra, additionally set subgraph size in bellman ford
void customize_bellman_ford_rebuild(crp::Graph *g, crp::OverlayStructure *overlay)
{
    auto init_algo = [&](int num_threads, uint32_t largest_cell) {
        return std::vector<BellmanFord>(num_threads, BellmanFord(largest_cell));
    };
    auto init_size = [&](BellmanFord &algo, NodeId subgraph_size, auto neighbors) {
        algo.set_number_of_nodes(subgraph_size);
    };
    auto one_to_all = [&](BellmanFord &algo, NodeId u, auto neighbors) { algo.compute_distance<false>(u, neighbors); };
    auto distance = [&](BellmanFord &algo, NodeId, NodeId to) { return algo.tentative_distance(to); };

    customize_shortest_path_rebuild(g, overlay, init_algo, init_size, one_to_all, distance);
}

void customize_floyd_warshall_rebuild(crp::Graph *g, crp::OverlayStructure *overlay)
{
    auto init_algo = [&](int num_threads, uint32_t largest_cell) {
        return std::vector<FloydWarshall>(num_threads, FloydWarshall(largest_cell));
    };
    auto init_size = [&](FloydWarshall &algo, NodeId subgraph_size, auto neighbors) {
        algo.set_number_of_nodes(subgraph_size);
        algo.compute_all_distances(neighbors);
    };
    auto one_to_all = [&](FloydWarshall, NodeId, auto) { return; };
    auto distance = [&](FloydWarshall &algo, NodeId u, NodeId v) { return algo.get_distance(u, v); };

    customize_shortest_path_rebuild(g, overlay, init_algo, init_size, one_to_all, distance);
}

void customize_bf_simd_rebuild(crp::Graph *g, crp::OverlayStructure *overlay)
{
    auto init_algo = [&](int num_threads, uint32_t largest_cell) {
        return std::vector<BellmanFordSIMD>(num_threads, BellmanFordSIMD(largest_cell));
    };
    auto compute_clique_entries = [&](LevelId level, CellId cellId, NodeId subgraph_size, auto neighbors,
                                      auto map_to_local, auto sp_algo) {
        const std::span<NodeId> border_nodes = overlay->get_border_nodes_for_cell(level, cellId);
        sp_algo.set_number_of_nodes(subgraph_size);
        for (auto it = border_nodes.begin(); it < border_nodes.end(); std::advance(it, SIMD_LEN))
        {
            std::array<NodeId, SIMD_LEN> start_nodes{};
            for (int i = 0; (i < SIMD_LEN) && (it + i) != border_nodes.end(); i++)
            {
                const NodeId u = *(it + i);
                start_nodes[i] = map_to_local(u);
            }
            sp_algo.compute_distance(start_nodes, neighbors);
            for (NodeId to : border_nodes)
            {
                const NodeId internal_to = map_to_local(to);
                const NodeId clique_to = overlay->get_internal_id(to, level);
                std::array<NodeId, SIMD_LEN> to_distances = sp_algo.tentative_distance(internal_to);
                for (int i = 0; (i < SIMD_LEN) && (it + i) != border_nodes.end(); i++)
                {
                    const NodeId u = *(it + i);
                    const NodeId clique_u = overlay->get_internal_id(u, level);
                    *overlay->get_distance(level, cellId, clique_u, clique_to) = to_distances[i];
                }
            }
        }
    };
    customize_rebuild_graph(g, overlay, init_algo, compute_clique_entries);
}

} // namespace crp
