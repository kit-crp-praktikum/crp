#include "crp/crp.h"
#include "data-types.h"

namespace crp
{

OverlayStructure::OverlayStructure(crp::Graph *g, RecursivePartition _partition) : partition(std::move(_partition))
{
    this->cliques.resize(partition.number_of_levels);

    const NodeId UNINITIALIZED_ID = g->num_nodes();

    this->node_id_on_level.resize(partition.number_of_levels, std::vector<NodeId>(g->num_nodes(), UNINITIALIZED_ID));
    this->border_nodes.resize(partition.number_of_levels);
    this->child_cell_ids.resize(partition.number_of_levels);

    int num_cells_in_level = 1;
    for (LevelId level = partition.number_of_levels - 1; level >= 0; level--)
    {
        num_cells_in_level *= partition.cells_per_level;

        cliques[level].resize(num_cells_in_level);
        border_nodes[level].resize(num_cells_in_level);
        child_cell_ids[level].resize(num_cells_in_level);

        const auto &try_add_node = [&](NodeId u) {
            if (node_id_on_level[level][u] == UNINITIALIZED_ID)
            {
                const auto cellU = partition.find_cell_for_node(u, level);
                node_id_on_level[level][u] = border_nodes[level][cellU].size();
                border_nodes[level][cellU].push_back(u);
            }
        };

        // Find border nodes for this level
        for (NodeId u = 0; u < (NodeId)g->num_nodes(); u++)
        {
            for (auto [v, _] : (*g)[u])
            {
                // Border nodes are nodes which have an edge to an adjacent cell.
                // To avoid duplicates, we only add undirected edges in one direction (u < v).
                if ((u < v || !g->get_edge(v, u)) && (partition.find_level_differing(u, v) >= level))
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

        // compute cellId of cells contained in level - 1
        // level 0 has no child cells -> empty
        if (level > 0)
        {
            const uint32_t bits_per_level = partition.get_bits_per_level();
            uint32_t bits = (partition.number_of_levels - level) * bits_per_level;
            CellId cells = partition.cells_per_level;
            CellId cells_in_level = num_cells_in_level;
            assert((cells & (cells - 1)) == 0);
            // if cells is a power of two -> all submasks are used
            for (CellId cell = 0; cell < cells_in_level; cell++)
            {
                for (CellId child_cell = 0; child_cell < cells; child_cell++)
                {
                    CellId c = (child_cell << bits) | cell;
                    child_cell_ids[level][cell].push_back(c);
                }
            }
        }
    }

    nodes_in_level_0.resize(num_cells_in_level);
    for (NodeId u = 0; u < (NodeId)g->num_nodes(); u++)
    {
        CellId cell = partition.find_cell_for_node(u, 0);
        nodes_in_level_0[cell].push_back(u);
    }
}

CellId OverlayStructure::get_cell_for_node(NodeId u, LevelId level)
{
    return partition.find_cell_for_node(u, level);
}

std::span<NodeId> OverlayStructure::get_border_nodes_for_cell(LevelId level, CellId cell)
{
    return border_nodes[level][cell];
}

std::span<CellId> OverlayStructure::get_child_cellIds(LevelId level, CellId cell)
{
    return child_cell_ids[level][cell];
}

std::span<NodeId> OverlayStructure::get_nodes_level0(CellId cell)
{
    return nodes_in_level_0[cell];
}

int OverlayStructure::num_cells_in_level(int level)
{
    return border_nodes[level].size();
}

NodeId OverlayStructure::get_internal_id(NodeId u, LevelId level)
{
    return node_id_on_level[level][u];
}

int OverlayStructure::get_number_of_levels()
{
    return partition.number_of_levels;
}

void OverlayStructure::remove_phantom_levels(int number_of_phantom_levels)
{
    partition.number_of_levels -= number_of_phantom_levels;

    // clear bitmask of phantom_levels in overlay partition
    for (uint32_t i = 0; i < partition.mask.size(); i++)
    {
        uint32_t bits = partition.get_bits_per_level() * partition.number_of_levels;
        uint32_t submask = ((uint32_t)1 << bits) - 1;
        partition.mask[i] = partition.mask[i] & submask;
    }

    // remove datastructures for phantom_levels
    cliques.erase(cliques.begin(), cliques.begin() + number_of_phantom_levels);
    node_id_on_level.erase(node_id_on_level.begin(), node_id_on_level.begin() + number_of_phantom_levels);
    border_nodes.erase(border_nodes.begin(), border_nodes.begin() + number_of_phantom_levels);

    // these datastructure are not needed for query
    child_cell_ids.clear();
    nodes_in_level_0.clear();
}

CRPAlgorithm::CRPAlgorithm(CRPAlgorithmParams params) : partition{params.number_of_levels, params.cells_per_level}
{
    this->params = params;
}

void OverlayStructure::precompute_cliquesT()
{
    this->cliquesT = this->cliques;
    for (LevelId level = 0; level < partition.number_of_levels; level++)
    {
        for (CellId cell = 0; cell < num_cells_in_level(level); cell++)
        {
            const NodeId n = cliques[level][cell].size();
            for (NodeId x = 0; x < n; x++)
            {
                for (NodeId y = 0; y < x; y++)
                {
                    std::swap(cliquesT[level][cell][x][y], cliquesT[level][cell][y][x]);
                }
            }
        }
    }
}

uint64_t OverlayStructure::total_border_nodes()
{
    uint64_t border_nodes = 0;
    for (LevelId level = 0; level < partition.number_of_levels; level++)
    {
        for (CellId cell = 0; cell < num_cells_in_level(level); cell++)
        {
            border_nodes += cliques[level][cell].size();
        }
    }
    return border_nodes;
}

uint64_t OverlayStructure::total_memory_bytes()
{
    uint64_t num_nodes = node_id_on_level[0].size();

    // node_id_on_level + nodes_in_level_0
    uint64_t memory_bytes = (partition.number_of_levels + 1) * num_nodes * sizeof(NodeId);

    for (LevelId level = 0; level < partition.number_of_levels; level++)
    {
        for (CellId cell = 0; cell < num_cells_in_level(level); cell++)
        {
            // clique + cliqueT + border_nodes
            memory_bytes +=
                (2 * sizeof(Distance) + sizeof(NodeId)) * cliques[level][cell].size() * cliques[level][cell].size();

            // child_cell_ids
            memory_bytes += sizeof(CellId) * child_cell_ids[level][cell].size();
        }
    }
    return memory_bytes;
}

uint64_t OverlayStructure::largest_cell()
{
    uint64_t largest_cell = nodes_in_level_0[0].size();
    for (CellId cell = 0; cell < num_cells_in_level(0); cell++)
    {
        largest_cell = std::max(largest_cell, nodes_in_level_0[cell].size());
    }
    return largest_cell;
}

} // namespace crp
