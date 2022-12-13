#include <iostream>
#include <vector>
#include <bitset>
#include "partitioner/bfs-partitioner.h"
#include "partitioner/rec-partitioner.h"
#include "graph.h"

int main2() {
    unsigned n = 32;
    std::vector<std::vector<std::pair<NodeId, Distance>>> gr(n);
    //path
    for(NodeId i = 0; i < n - 1; i++) {
        gr[i].push_back({i+1,1});
        gr[i + 1].push_back({i,1});
    }

    crp::Graph g(gr);

    std::vector<Distance> tmp_lat;
    partitioner::BfsPartitioner bfs;
    partitioner::RecPartitioner recPart(&bfs, &g);

    std::vector<NodeId> masks = recPart.partition_rec(4,2);

    for (unsigned i = 0; i < masks.size(); i++) {
        std::cout << "Node " << i << ": " <<  std::bitset<8>(masks[i]) << std::endl;
    }





    

}