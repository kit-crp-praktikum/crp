#include "graph.h"
#include "partitioner/bfs-partitioner.h"
#include "partitioner/rec-partitioner.h"
#include <bitset>
#include <iostream>
#include <vector>

int main()
{
    unsigned n = 32;
    std::vector<std::vector<std::pair<NodeId, Distance>>> gr(n);
    // path
    for (NodeId i = 0; i < n - 1; i++)
    {
        gr[i].push_back({i + 1, 1});
        gr[i + 1].push_back({i, 1});
    }

    partitioner::BfsPartitioner bfs;
    partitioner::RecPartitioner recPart(&bfs);

    std::vector<NodeId> *masks = recPart.partition_rec(3, 2, &gr);

    for (unsigned i = 0; i < masks->size(); i++)
    {
        std::cout << "Node " << i << ": " << (std::bitset<8>)((*masks)[i]) << std::endl;
    }
}