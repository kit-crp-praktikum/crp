#include "kahip-wrapper.hpp"
#include "crp/crp.h"
#include "data-types.h"
#include "graph.h"
#include "partitioner/preprocessing.hpp"
#include "subgraph.h"
#include <numeric>

#if !ENABLE_KAHIP
crp::RecursivePartitionMask partitioner::kahip_partition_graph(crp::Graph *g, int nr_levels, int nr_cells_per_level)
{
    std::cerr << "CRP built without KaHIP support!" << std::endl;
    std::exit(-1);
}

#else
#include "kaHIP_interface.h"
static void recursive_partition(partitioner::Subgraph graph, crp::RecursivePartitionMask &out_mask, int nr_levels,
                                int nr_cells_per_level)
{
    if (nr_levels == 0)
    {
        return;
    }

    crp::Graph g{graph.graph};

    int n = g.num_nodes();
    double imbalance = 0.05;
    int edgecut = 0; // idk what this is about

    std::vector<int> partition(n);

    kaffpa(&n, NULL, (int *)g.first_out.data(), NULL, (int *)g.head.data(), &nr_cells_per_level, &imbalance, false, 0,
           ECO, &edgecut, partition.data());

    int bits_per_level = 32 - __builtin_clz(nr_cells_per_level - 1);

    auto subs = partitioner::generate_subgraphs(graph, partition, nr_cells_per_level);
#pragma omp parallel for
    for (int i = 0; i < nr_cells_per_level; i++)
    {
        recursive_partition(subs[i], out_mask, nr_levels - 1, nr_cells_per_level);
        for (NodeId x = 0; x < subs[i].graph.size(); x++)
        {
            auto &mask = out_mask[subs[i].mapping[x]];
            mask <<= bits_per_level;
            mask |= i;
        }
    }
}

crp::RecursivePartitionMask partitioner::kahip_partition_graph(crp::Graph *g, int nr_levels, int nr_cells_per_level)
{
    std::vector<NodeId> mapping(g->num_nodes());
    std::iota(mapping.begin(), mapping.end(), 0u);

    partitioner::Subgraph graph{
        .graph = make_undirected(*g),
        .mapping = mapping,
        .geo_data = {},
    };

    crp::RecursivePartitionMask mask(g->num_nodes());
    recursive_partition(graph, mask, nr_levels, nr_cells_per_level);
    return mask;
}
#endif
