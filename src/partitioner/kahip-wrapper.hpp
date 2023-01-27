#pragma once

#include "crp/crp.h"
#include "graph.h"

namespace partitioner
{
crp::RecursivePartitionMask kahip_partition_graph(crp::Graph *g, int nr_levels, int nr_cells_per_level);
}
