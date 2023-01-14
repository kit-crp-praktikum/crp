#include <bitset>
#include <iostream>
#include <vector>

#include "graph.h"
#include "partitioner/bfs-partitioner.h"
#include "partitioner/rec-partitioner.h"
#include "tests/grid-graph.hpp"


int main()
{
    uint32_t n = 8;
    crp::AdjacencyList gr = generate_grid_graph(n);
    partitioner::GeoData geo_data = generate_grid_graph_embedding(n);

    int cells_per_level = 3;
    int number_of_levels = 2;
    partitioner::BfsPartitioner bfs;
    partitioner::RecPartitioner recPart(bfs, cells_per_level, number_of_levels);

    std::vector<NodeId> masks = recPart.partition_rec(gr, geo_data);
    crp::RecursivePartition partition = {number_of_levels, cells_per_level, masks};
    std::cout << "NodeIds \n";
    dump_grid_graph_node_ids(n);
    std::cout << "Partition \n";
    grid_graph_print_recursive_partition(n, partition );
}