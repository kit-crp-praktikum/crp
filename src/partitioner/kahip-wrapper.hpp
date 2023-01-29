#pragma once

#include "crp/crp.h"
#include "graph.h"

namespace partitioner
{
enum class KaHIPMode
{
    FAST = 0,
    ECO = 1,
    STRONG = 2,
    FASTSOCIAL = 3,
    ECOSOCIAL = 4,
    STRONGSOCIAL = 5,
};

struct KaHIPParameters
{
    double imbalance;
    KaHIPMode mode;
    int nr_levels;
    int nr_cells_per_level;
};

crp::RecursivePartitionMask kahip_partition_graph(crp::Graph *g, const KaHIPParameters &params);
} // namespace partitioner
