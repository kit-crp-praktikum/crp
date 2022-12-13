#pragma once

#include <vector>
#include "src/data-types.h"
#include "geo-data.h"
#include "src/graph.h"

/**
 * Partition graph by running a BFS until the half of the nodes are reached.
*/

namespace partitioner
{
class BfsPartitioner {
    public:
    BfsPartitioner();

    /**
     * Returns partition as bool vector of length n
     * Starts BFS from highest latitude point if info available
    */
    std::vector<bool> partition (crp::Graph* g, std::vector<bool>* cont_nodes);

    NodeId findStart(std::vector<bool>* visited, std::vector<bool>* cont_nodes);

    public:
    GeoData gd;

};
}