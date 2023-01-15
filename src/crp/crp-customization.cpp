#include "algorithms/dijkstra.hpp"
#include "crp.h"

#include <iostream>

namespace crp
{
void CRPAlgorithm::customize()
{
    this->reverse = g->reversed();
    this->params.customizer(g, overlay.get());
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
                const CellId upper_cell_u = overlay->get_cell_for_node(u, level);
                // Iterate graph edges
                for (auto [to, weight] : (*g)[u])
                {
                    if (upper_cell_u == overlay->get_cell_for_node(to, level))
                    {
                        f(to, weight);
                    }
                }

                // Iterate clique
                const LevelId prev_level = level - 1;
                const CellId lower_cellId = overlay->get_cell_for_node(u, prev_level);
                const NodeId internal_u = overlay->get_internal_id(u, prev_level);
                const std::span<NodeId> neighbors = overlay->get_border_nodes_for_cell(prev_level, lower_cellId);
                // invalid ID -> is no border node
                if (internal_u == g->num_nodes())
                {
                    return;
                }
                for (NodeId internal_to = 0; internal_to < neighbors.size(); internal_to++)
                {
                    f(neighbors[internal_to],
                      *overlay->get_distance(prev_level, lower_cellId, internal_u, internal_to));
                }
            };
            compute_clique(level, cellId, neighbors_in_cell);
        }
    }
}

void customize_with_dijkstra(crp::Graph *g, crp::OverlayStructure *overlay)
{
    Dijkstra one_to_all(g->num_nodes());
    auto compute_clique = [&](LevelId level, CellId cellId, auto neighbors_in_cell)
    { 
        const std::span<NodeId> border_nodes = overlay->get_border_nodes_for_cell(level, cellId);
        for (NodeId u : border_nodes)
        {
            one_to_all.compute_distance<false>(u, neighbors_in_cell);
            const NodeId internal_u = overlay->get_internal_id(u, level);
            for (NodeId to : border_nodes)
            {
                const NodeId internal_to = overlay->get_internal_id(to, level);
                *overlay->get_distance(level, cellId, internal_u, internal_to) = one_to_all.tentative_distance(to);
            }
        }
    };
    generic_customize(g, overlay, compute_clique);
}

} // namespace crp
