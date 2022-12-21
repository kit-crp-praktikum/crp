#include "crp/crp.h"
#include "data-types.h"

namespace crp {

OverlayStructure::OverlayStructure(crp::Graph* g, RecursivePartition _partition)
{
    this->partition = std::move(_partition);

    this->cliques.resize(partition.number_of_levels);

    const NodeId UNINITIALIZED_ID = g->num_nodes();

    this->node_id_on_level.resize(partition.number_of_levels,
        std::vector<NodeId>(g->num_nodes(), UNINITIALIZED_ID));
    this->border_nodes.resize(partition.number_of_levels);

    int num_cells_in_level = 1;
    for (LevelId level = 0; level < partition.number_of_levels; level++)
    {
        num_cells_in_level *= partition.cells_per_level;

        cliques[level].resize(num_cells_in_level);
        border_nodes[level].resize(num_cells_in_level);

        const auto& try_add_node = [&] (NodeId u)
        {
            if (node_id_on_level[level][u] == UNINITIALIZED_ID)
            {
                const auto cellU = partition.find_cell_for_node(u, level);
                node_id_on_level[level][u] = border_nodes[level][cellU].size();
                border_nodes[level][cellU].push_back(cellU);
            }
        };

        // Find border nodes for this level
        for (NodeId u = 0; u < g->num_nodes(); u++)
        {
            for (auto [v, _] : (*g)[u])
            {
                // Border nodes are nodes which have an edge to an adjacent cell.
                // To avoid duplicates, we only add undirected edges in one direction (u < v).
                if ((u < v || !g->get_edge(v, u)) &&
                    (partition.find_level_differing(u, v) == level))
                {
                    try_add_node(u);
                    try_add_node(v);
                }
            }
        }

        for (int cell = 0; cell < num_cells_in_level; cell++)
        {
            int nr_nodes = border_nodes[level][cell].size();
            cliques[level][cell] = Clique(nr_nodes, std::vector<Distance>(nr_nodes, INF));
        }
    }
}

CellId OverlayStructure::get_cell_for_node(NodeId u, LevelId level)
{
    return partition.find_cell_for_node(u, level);
}

std::span<NodeId> OverlayStructure::get_border_nodes_for_cell(
    LevelId level, CellId cell)
{
    return border_nodes[level][cell];
}

Distance* OverlayStructure::get_distance(
    LevelId level, CellId cell, NodeId a, NodeId b)
{
    return &cliques[level][cell][a][b];
}

int OverlayStructure::num_cells_in_level(int level)
{
    return border_nodes[level].size();
}

NodeId OverlayStructure::get_internal_id(NodeId u, LevelId level)
{
    return node_id_on_level[level][u];
}

} // namespace crp
