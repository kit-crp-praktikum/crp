#pragma once

#include <algorithm>
#include <vector>
#include <iostream>

#include "graph.h"


namespace crp
{

class PartitionStats {
    public:
    PartitionStats(crp::Graph* g,  crp::OverlayStructure* overlay);

    float get_avg_border_nodes(int level);

    float get_avg_cell_size(int level);

    // returns ratio of quartiles: size(0.75)/size(0.25)
    float get_cell_size_imbalance(int level);

    void print_stats();

    private:
    crp::Graph *g;
    crp::OverlayStructure *overlay;
};

}