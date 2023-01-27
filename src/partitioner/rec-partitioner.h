#pragma once

#include <algorithm>
#include <numeric>
#include <vector>

#include "bipartitioner.h"
#include "src/data-types.h"
#include "subgraph.h"

/**
 * Recursive partitioning using 2-partitioner interface
 * Bit mask for each node is created
 * Lowest lewel bits are being pushed to the end of the mask
 */

using ClusterId = uint32_t;

namespace partitioner
{
class RecPartitioner
{
  public:
    RecPartitioner(BiPartitioner &_bi_partitioner, uint32_t _cells_per_level, uint32_t _number_of_levels)
        : bi_partitioner(_bi_partitioner), cells_per_level(_cells_per_level), number_of_levels(_number_of_levels)
    {
    }

    /**
     * Partition graph recursively.
     * Creates bit masks for partitions.
     */
    std::vector<NodeId> partition_rec(const crp::AdjacencyList &graph, const GeoData &geo_data)
    {
        std::vector<ClusterId> masks(graph.size(), 0);
        std::vector<NodeId> map_default(graph.size());
        std::iota(map_default.begin(), map_default.end(), 0);

        // copy original graph, since we don't want to change input graph.

        Subgraph original_graph{
            .graph = graph,
            .mapping = std::move(map_default),
            .geo_data = geo_data,
        };

        partition_rec_impl(original_graph, masks, number_of_levels);
        return masks;
    }

  private:
    /**
     * Returns a bitmask for each node indicating which cells it belongs to.
     */
    void partition_rec_impl(Subgraph &graph, std::vector<ClusterId> &masks, uint32_t current_level)
    {
        if (current_level == 0)
            return;

        // at the beginning we only have one cell with one mapping to the original graph
        std::vector<Subgraph> subgraphs;
        subgraphs.push_back(std::move(graph));

        // divide largest cell in two until there are cells_per_level cells
        while (subgraphs.size() < cells_per_level)
        {
            // partition largest cell
            uint32_t divide_pos = find_largest_cell(subgraphs);

            // swap largest cell to last position
            std::swap(subgraphs[divide_pos], subgraphs.back());
            std::vector<bool> bipartition = bi_partitioner.partition(subgraphs.back().graph, subgraphs.back().geo_data);

            // calculate new cells and mappings for last graph in the list and delete it
            auto new_cells = generate_subgraphs(subgraphs.back(), bipartition, 2);
            subgraphs.pop_back();
            subgraphs.insert(subgraphs.end(), new_cells.begin(), new_cells.end());
        }

        update_masks(masks, subgraphs, current_level);

        // divide cells into cells on lower level recursively
        for (uint32_t i = 0; i < cells_per_level; i++)
        {
            partition_rec_impl(subgraphs[i], masks, current_level - 1);
        }
        return;
    }

    uint32_t find_largest_cell(std::vector<Subgraph> &subgraphs)
    {
        uint32_t divide_pos = 0;
        for (uint32_t i = 1; i < subgraphs.size(); i++)
        {
            if (subgraphs[i].graph.size() > subgraphs[divide_pos].graph.size())
                divide_pos = i;
        }
        return divide_pos;
    }

    void update_masks(std::vector<ClusterId> &masks, std::vector<Subgraph> &subgraphs, uint32_t level)
    {
        const uint32_t bits_per_level = 32 - __builtin_clz(cells_per_level - 1);
        const uint32_t shift = bits_per_level * (number_of_levels - level);
        // update masks for current layer's cells
        for (ClusterId part_id = 0; part_id < subgraphs.size(); part_id++)
        {
            for (auto v_orig : subgraphs[part_id].mapping)
            {
                // v_orig: original index in graph
                masks[v_orig] = masks[v_orig] | (part_id << shift);
            }
        }
    }

  private:
    BiPartitioner &bi_partitioner;
    uint32_t cells_per_level;
    uint32_t number_of_levels;
};
} // namespace partitioner
