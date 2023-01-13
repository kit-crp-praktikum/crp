#pragma once

#include "bfs-partitioner.h"
#include "src/data-types.h"
#include <vector>

/**
 * Recursive partitioning using 2-partitioner interface
 * Bit mask for each node is created
 * Lowest lewel bits are being pushed to the end of the mask
 */

using ClusterId = uint32_t;

namespace partitioner
{
template <class BiPartitioner> class RecPartitioner
{
  public:
    RecPartitioner(BiPartitioner *bip) : biPartitioner(bip)
    {
    }

    /**
     * Partition graph recursively
     * Create bit masks for partitions
     * @param cells_per_level: # of cells on each level
     * @param levels: # of levels
     * @param g: graph to be partitioned
     * @return vector of bit masks for each node
     */
    std::vector<NodeId> *partition_rec(int cells_per_level, int levels,
                                       std::vector<std::vector<std::pair<NodeId, Distance>>> *g)
    {
        masks.clear();
        masks.resize(g->size());

        // shift level info into masks
        bits_for_mask = 0;
        int shifted = cells_per_level;
        while (shifted > 0)
        {
            shifted = (shifted >> 1);
            bits_for_mask++;
        }

        std::vector<NodeId> map_default(g->size());
        for (NodeId i = 0; i < g->size(); i++)
        {
            map_default[i] = i;
        }

        partition_rec_impl(cells_per_level, levels, g, &map_default);
        return &masks;
    }

  private:
    /**
     * Returns a bitmask for each node indicating which cells it belongs to.
     * Private version for recursive calls
     * @param cells_per_level: # of cells on each level
     * @param levels: # of levels
     * @param g: graph to be partitioned
     */
    void partition_rec_impl(int cells_per_level, int levels, std::vector<std::vector<std::pair<NodeId, Distance>>> *g,
        std::vector<NodeId> * map_default)
    {
        if (levels == 0)
            return;

        // variables for graphs and mappings of the new partitiones
        std::vector<NodeId> map_new_orig0[cells_per_level], map_new_orig1[cells_per_level];
        std::vector<std::vector<std::pair<NodeId, Distance>>> cell0[cells_per_level], cell1[cells_per_level];

        // contains each current cell's graph on this level, and node id mapping to original id-s
        std::vector<std::vector<std::vector<std::pair<NodeId, Distance>>> *> cells;
        std::vector<std::vector<NodeId> *> mappings;

        mappings.push_back(map_default);
        cells.push_back(g);

        // divide largest cell in two until there are cells_per_level cells
        for (int cell_cnt = 1; cell_cnt < cells_per_level; cell_cnt++)
        {
            // find largest cell
            unsigned divide_pos = 0;
            for (unsigned i = 1; i < cells.size(); i++)
            {
                if (cells[i]->size() > cells[divide_pos]->size())
                    divide_pos = i;
            }
            std::vector<std::vector<std::pair<NodeId, Distance>>> *cell_to_divide = cells[divide_pos];

            // partition largest cell
            std::vector<bool>* nodes01 = (*biPartitioner).partition(cell_to_divide);

            // calculate new cells and mappings
            divide_graph(mappings[divide_pos], nodes01, cell_to_divide, &map_new_orig0[cell_cnt],
                         &map_new_orig1[cell_cnt], &cell0[cell_cnt], &cell1[cell_cnt]);

            // remove old cell, add two new cells
            cells.erase(cells.begin() + divide_pos);
            mappings.erase(mappings.begin() + divide_pos);
            cells.push_back(&cell0[cell_cnt]);
            mappings.push_back(&map_new_orig0[cell_cnt]);
            cells.push_back(&cell1[cell_cnt]);
            mappings.push_back(&map_new_orig1[cell_cnt]);
        }

        // update masks for current layer's cells
        for (ClusterId part_id = 0; part_id < mappings.size(); part_id++)
        {
            std::vector<NodeId> *part_mapping = mappings[part_id];
            for (auto v_orig : (*part_mapping))
            {   
                // v_orig: original index in graph
                masks[v_orig] = masks[v_orig] << bits_for_mask;
                masks[v_orig] = masks[v_orig] | part_id;
            }
            
        }

        // divide cells into cells on lower level recursively
        for (int i = 0; i < cells_per_level; i++)
        {
            partition_rec_impl(cells_per_level, levels - 1, cells[i], mappings[i]);
        }
        return;
    }

    /**
     * divide current graph into two smaller ones
     * there are three kinds of graphs, each with its own node labeling:
     * - original: the graph at the very beginning
     * - actual: graph which is being divided
     * - new: the newly created two graphs
     * @param map_act_orig: maps actual graph indices to original ones
     * @param nodes01: nodes of two partitions coded as 0-1 vector
     * @param g: graph to divide (actual graph)
     * @return multiple values by pointers: map_new_orig_0, map_new_orig_1, g0, g1
     *
     */
    void divide_graph(std::vector<NodeId> *map_act_orig, std::vector<bool> *part_nodes,
                      std::vector<std::vector<std::pair<NodeId, Distance>>> *g, std::vector<NodeId> *map_new_orig_0,
                      std::vector<NodeId> *map_new_orig_1, std::vector<std::vector<std::pair<NodeId, Distance>>> *g0,
                      std::vector<std::vector<std::pair<NodeId, Distance>>> *g1)
    {
        // create new_orig mapping
        std::vector<NodeId> map_act_new_0(g->size());
        std::vector<NodeId> map_act_new_1(g->size());

        for (unsigned i = 0; i < g->size(); i++)
        {
            map_act_new_0[i] = invalid_id;
            map_act_new_1[i] = invalid_id;
        }

        unsigned g0_size = std::count(part_nodes->begin(), part_nodes->end(), false);
        unsigned g1_size = part_nodes->size() - g0_size;

        (*map_new_orig_0).resize(g0_size);
        (*map_new_orig_1).resize(g1_size);

        NodeId new_node_id0 = 0;
        NodeId new_node_id1 = 0;

        // create new->orig mappings
        for (NodeId node = 0; node < part_nodes->size(); node++)
        {
            if (((*part_nodes)[node] == false))
            {
                // node in partition 0
                (*map_new_orig_0)[new_node_id0] = (*map_act_orig)[node];
                map_act_new_0[node] = new_node_id0;
                new_node_id0++;
            }
            else
            {
                // node in partition 1
                (*map_new_orig_1)[new_node_id1] = (*map_act_orig)[node];
                map_act_new_1[node] = new_node_id1;
                new_node_id1++;
            }
        }

        // create new graphs
        g0->resize(g0_size);
        g1->resize(g1_size);
        // find u->v edge in actual graph, add it to new graph
        for (NodeId u_act = 0; u_act < part_nodes->size(); u_act++) {
            if ((*part_nodes)[u_act] == false) {
                // u_act in g0
                NodeId u_new = map_act_new_0[u_act];
                for (auto [v_act, weight] : (*g)[u_act]) {
                    NodeId v_new = map_act_new_0[v_act];
                    // check if v is also in g0's partition
                    if (v_new != invalid_id)
                    {
                        (*g0)[u_new].push_back({v_new, weight});
                    }
                }
            }
            else {
                // u_act in g1
                NodeId u_new = map_act_new_1[u_act];
                for (auto [v_act, weight] : (*g)[u_act]) {
                    NodeId v_new = map_act_new_1[v_act];
                    // check if v is also in g1's partition
                    if (v_new != invalid_id)
                    {
                        (*g1)[u_new].push_back({v_new, weight});
                    }
                }
            }
        }
    }

  private:
    BiPartitioner *biPartitioner;
    std::vector<ClusterId> masks;
    int bits_for_mask;
};
} // namespace partitioner