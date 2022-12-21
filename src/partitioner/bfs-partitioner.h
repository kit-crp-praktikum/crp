#pragma once

#include "geo-data.h"
#include "src/data-types.h"
#include <unordered_set>
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
     * Returns two partitions as pair of node sets
     * Starts BFS from highest latitude point if info available
     */
    std::pair<std::unordered_set<NodeId> *, std::unordered_set<NodeId> *> partition(
        std::vector<std::vector<std::pair<NodeId, Distance>>> *g);

  private:
    /**
     * Find new start node if component finished
     */
    NodeId findStart(std::vector<bool> *visited);

  private:
    GeoData gd;
    NodeId lastStart = 0;
    std::unordered_set<NodeId> result0, result1;
};
} // namespace partitioner