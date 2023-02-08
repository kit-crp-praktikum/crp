#pragma once

#include <algorithm>
#include <vector>
#include <iostream>

#include "crp.h"

#include "partition-stats.h"

namespace crp
{

PartitionStats::PartitionStats(crp::Graph* _g,  crp::OverlayStructure* _overlay) : g(_g), overlay(_overlay) 
    {
    }

    float PartitionStats::get_avg_border_nodes(int level) 
    {
        int num_cells = overlay->num_cells_in_level(level);
        int num_border_nodes = 0;
        for (CellId cell = 0; cell < num_cells; cell++) {
            num_border_nodes += overlay->get_border_nodes_for_cell(level, cell).size();
        }
        return (float)num_border_nodes / num_cells;
    }

    float PartitionStats::get_avg_cell_size(int level) 
    {
        int num_nodes = g->num_nodes();
        int num_cells = overlay->num_cells_in_level(level);
        return (float)num_nodes / num_cells;
    }

    // returns ratio of quartiles: size(0.75)/size(0.25)
    float PartitionStats::get_cell_size_imbalance(int level) 
    {
        int num_cells = overlay->num_cells_in_level(level);
        std::vector<int> cell_sizes(num_cells);
        for (CellId i = 0; i < num_cells; i++) {
            cell_sizes[i] = overlay->get_border_nodes_for_cell(level, i).size();
        }
        std::sort(cell_sizes.begin(), cell_sizes.end());
        return (float)cell_sizes[num_cells * 3/4] / cell_sizes[num_cells /4];
    }

    void PartitionStats::print_stats()
    {
        std::cout << "avg border nodes:" << get_avg_border_nodes(0);
        std::cout << "\navg cell size: " << get_avg_cell_size(0);
        std::cout << "\nimbalance: " << get_cell_size_imbalance(0) << "\n";
    }

}