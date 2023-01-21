#include "algorithms/bellman_ford.hpp"
#include "algorithms/dijkstra.hpp"
#include "algorithms/floyd_warshall.hpp"
#include "crp.h"

#include <iostream>

namespace crp
{
void CRPAlgorithm::customize()
{
    this->reverse = g->reversed();
    // This is the overlay on the recursive partition with phantom_levels.
    this->overlay = std::make_unique<OverlayStructure>(this->g, this->partition);
    this->params.customizer(g, overlay.get());    
    this->overlay->remove_phantom_levels(this->params.number_of_phantom_levels);
}

void generic_customize(crp::Graph *g, crp::OverlayStructure *overlay, auto compute_clique)
{
    // compute level 0 cliques from original graph
    LevelId level = 0;
    for (CellId cellId = 0; cellId < overlay->num_cells_in_level(level); cellId++)
    {
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
        compute_clique(level, cellId, neighbors_in_cell);
    }

    // compute level i cliques from level i - 1 cliques
    for (LevelId level = 1; level < overlay->get_number_of_levels(); level++)
    {   
        for (CellId cellId = 0; cellId < overlay->num_cells_in_level(level); cellId++)
        {
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
            compute_clique(level, cellId, neighbors_in_cell);
        }
    }
}


void customize_with_shortest_path(crp::Graph *g, crp::OverlayStructure *overlay, auto init_algo,
                                     auto process_node, auto distance)
{
    auto algo = init_algo();

    auto compute_clique = [&](LevelId level, CellId cellId, auto neighbors) {
        const std::span<NodeId> border_nodes = overlay->get_border_nodes_for_cell(level, cellId);
        for (NodeId u : border_nodes)
        {
            process_node(algo, u, neighbors, level, cellId);
            const NodeId internal_u = overlay->get_internal_id(u, level);
            for (NodeId to : border_nodes)
            {
                const NodeId internal_to = overlay->get_internal_id(to, level);
                *overlay->get_distance(level, cellId, internal_u, internal_to) = distance(algo, to);
            }
        }
    };
    generic_customize(g, overlay, compute_clique);
}

void customize_with_dijkstra(crp::Graph *g, crp::OverlayStructure *overlay)
{
    auto init_algo = [&]() { return Dijkstra(g->num_nodes()); };
    auto process_node = [&](Dijkstra &algo, NodeId u, auto neighbors, LevelId, CellId) { algo.compute_distance<false>(u, neighbors); };
    auto distance = [&](Dijkstra &algo, NodeId to) { return algo.tentative_distance(to); };

    customize_with_shortest_path(g, overlay, init_algo, process_node, distance);
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
    // we don't set bellman ford size to subgraph size, but it will probably converge earlier -> worst case 1 iteration more than necessary 
    auto init_algo = [&]() { return BellmanFord(g->num_nodes()); };
    auto process_node = [&](BellmanFord &algo, NodeId u, auto neighbors, LevelId level, CellId cellId) 
    {   
        auto nodes = [&](auto f) {for_all_nodes(level, cellId, f);};
        algo.generic_compute_distance<false>(u, neighbors, nodes); 
    };
    auto distance = [&](BellmanFord &algo, NodeId to) { return algo.tentative_distance(to); };
    customize_with_shortest_path(g, overlay, init_algo, process_node, distance);
}

void customize_rebuild_graph(crp::Graph *g, crp::OverlayStructure *overlay, auto compute_clique_entries)
{
    crp::Graph cell_graph;
    std::vector<NodeId> original_to_local_node(g->num_nodes());

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

    auto build_graph = [&](LevelId level, CellId cellId, auto neighbors_in_cell) {
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

    auto compute_clique = [&](LevelId level, CellId cellId, auto neighbors_in_cell) {
        build_graph(level, cellId, neighbors_in_cell);
        auto neighbor_in_cell_graph = [&](NodeId v, auto f) {
            for (auto [u, weight] : cell_graph[v])
            {
                f(u, weight);
            }
        };
        auto map_to_local = [&](NodeId v) { return original_to_local_node[v]; };
        compute_clique_entries(level, cellId, cell_graph.num_nodes(), neighbor_in_cell_graph, map_to_local);
    };
    generic_customize(g, overlay, compute_clique);
}

void customize_shortest_path_rebuild(crp::Graph *g, crp::OverlayStructure *overlay, auto init_algo, auto init_size,
                                     auto process_node, auto distance)
{
    auto algo = init_algo();

    auto compute_clique_entries = [&](LevelId level, CellId cellId, NodeId subgraph_size, auto neighbors,
                                      auto map_to_local) {
        init_size(algo, subgraph_size, neighbors);
        const std::span<NodeId> border_nodes = overlay->get_border_nodes_for_cell(level, cellId);
        for (NodeId u : border_nodes)
        {
            const NodeId internal_u = map_to_local(u);
            process_node(algo, internal_u, neighbors);
            for (NodeId to : border_nodes)
            {
                const NodeId internal_to = map_to_local(to);
                const NodeId clique_u = overlay->get_internal_id(u, level);
                const NodeId clique_to = overlay->get_internal_id(to, level);
                *overlay->get_distance(level, cellId, clique_u, clique_to) = distance(algo, internal_u, internal_to);
            }
        }
    };
    customize_rebuild_graph(g, overlay, compute_clique_entries);
}

void customize_dijkstra_rebuild(crp::Graph *g, crp::OverlayStructure *overlay)
{

    auto init_algo = [&]() { return Dijkstra(g->num_nodes()); };
    auto init_size = [&](Dijkstra, NodeId, auto) { return; };
    auto process_node = [&](Dijkstra &algo, NodeId u, auto neighbors) { algo.compute_distance<false>(u, neighbors); };
    auto distance = [&](Dijkstra &algo, NodeId, NodeId to) { return algo.tentative_distance(to); };

    customize_shortest_path_rebuild(g, overlay, init_algo, init_size, process_node, distance);
}

// same code as dijkstra, additionally set subgraph size in bellman ford
void customize_bellman_ford_rebuild(crp::Graph *g, crp::OverlayStructure *overlay)
{
    auto init_algo = [&]() { return BellmanFord(g->num_nodes()); };
    auto init_size = [&](BellmanFord &algo, NodeId subgraph_size, auto neighbors) {
        algo.set_number_of_nodes(subgraph_size);
    };
    auto process_node = [&](BellmanFord &algo, NodeId u, auto neighbors) {
        algo.compute_distance<false>(u, neighbors);
    };
    auto distance = [&](BellmanFord &algo, NodeId, NodeId to) { return algo.tentative_distance(to); };

    customize_shortest_path_rebuild(g, overlay, init_algo, init_size, process_node, distance);
}

void customize_floyd_warshall_rebuild(crp::Graph *g, crp::OverlayStructure *overlay)
{
    auto init_algo = [&]() { return FloydWarshall(FLOYD_WARSHALL_MAX_N); };
    auto init_size = [&](FloydWarshall &algo, NodeId subgraph_size, auto neighbors) {
        algo.set_number_of_nodes(subgraph_size);
        algo.compute_all_distances(neighbors);
    };
    auto process_node = [&](FloydWarshall, NodeId, auto) { return; };
    auto distance = [&](FloydWarshall &algo, NodeId u, NodeId v) { return algo.get_distance(u, v); };

    customize_shortest_path_rebuild(g, overlay, init_algo, init_size, process_node, distance);
}

} // namespace crp
