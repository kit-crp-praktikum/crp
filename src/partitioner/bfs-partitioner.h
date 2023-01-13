#pragma once

#include "geo-data.h"
#include "src/data-types.h"

#include <vector>

/**
 * Partition graph by running a BFS until the half of the nodes are reached.
 */

namespace partitioner
{
class BfsPartitioner
{
  public:
    BfsPartitioner();

    /**
     * Returns two partitions as bool vector
     * Starts BFS from highest latitude point if info available
     */
    std::vector<bool>* partition(
        std::vector<std::vector<std::pair<NodeId, Distance>>> *g);

  private:
    /**
     * Find new start node if component finished
     */
    NodeId findStart(std::vector<bool> *visited);

  private:
    GeoData gd;
    NodeId lastStart = 0;
    std::vector<bool> visited;
};
} // namespace partitioner