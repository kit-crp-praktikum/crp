#pragma once

#include <vector>
#include "src/data-types.h"
#include "src/graph.h"
#include "bfs-partitioner.h"

/**
 * Recursive partitioning using 2-partitioner interface
 * Bit mask for each node is created 
 * Lowest lewel bits are being pushed to the end of the mask
*/

namespace partitioner
{
template<class BiPartitioner>
class RecPartitioner {
    public:
    RecPartitioner(BiPartitioner* bip, crp::Graph* _g) : biPartitioner(bip), g(_g) 
    {
        masks.resize(g->num_nodes());
    }


    /**
     * Returns a bitmask for each node indicating which cells it belongs to.
     * @param cells_per_level: # of cells on each level, power of two!
     * @param levels: # of levels
     * @param cont_nodes: indicatir for node set of graph which has to be partitioned in two
    */
    std::vector<NodeId> partition_rec(int cells_per_level, int levels)
    {
        std::vector<bool> cont_nodes(g->num_nodes(), 1);
        partition_rec(cells_per_level, cells_per_level, levels, &cont_nodes);
        return masks;
    }

    private: 
    /**
     * Private version for recursive calls
     * @param cells_per_level: # of cells on each level
     * @param cells_remaining: # of cells on this level must be subdivided into this number of cells
    */
    void partition_rec(int cells_per_level, int cells_remaining, int levels, std::vector<bool>* cont_nodes)
    {
        if (levels == 0) return;

        // cont_part0 and cont_part1 are the cont_nodes of the two new partitions
        std::vector<bool> cont_part0 = (*biPartitioner).partition(g, cont_nodes);
        std::vector<bool> cont_part1 = cont_part0;
        cont_part1.flip();
        for (unsigned i = 0; i < cont_part1.size(); i++) 
        {
            cont_part1[i] = cont_part1[i] & (*cont_nodes)[i];
        }

        // add 0 to end of cont_part0 masks
        // add 1 to end of cont_part1 masks
        for (int i = 0; i < g->num_nodes(); i++)
        {
            if (cont_part0[i]) {
                masks[i] = masks[i]<<1;
            }
            else if (cont_part1[i]) {
                masks[i] = masks[i] << 1;
                masks[i] = masks[i] | 1;
            }
        }

        if (cells_remaining > 2)
        {
            // stay on level, current cells must be further divided
            partition_rec(cells_per_level, cells_remaining/2, levels, &cont_part0);
            partition_rec(cells_per_level, cells_remaining/2, levels, &cont_part1);
        }
        else 
        {
            // this level is done, step to next level
            partition_rec(cells_per_level, cells_per_level, levels-1, &cont_part0);
            partition_rec(cells_per_level, cells_per_level, levels-1, &cont_part1);
        }

        return;
    }

    
    private:
    BiPartitioner* biPartitioner;
    crp::Graph* g;
    std::vector<NodeId> masks;

};
}