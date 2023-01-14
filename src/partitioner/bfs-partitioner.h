#pragma once

#include <vector>
#include <algorithm>
#include <queue>
#include <iostream>

#include "geo-data.h"
#include "src/data-types.h"
#include "graph.h"
#include "geo-data.h"
#include "bipartitioner.h"


/**
 * Partition graph by running a BFS until the half of the nodes are reached.
 */
namespace partitioner
{
class BfsPartitioner : public BiPartitioner
{
  public:
    BfsPartitioner();

    /**
     * Returns two partitions as bool vector
     * Starts BFS from highest latitude point if info available
     */
    std::vector<bool> partition(crp::AdjacencyList &graph, GeoData &geo_data);

  private:
    /**
     * Find new start node if component finished
     */
    NodeId findStart(std::vector<bool> &visited, GeoData &geo_data);

  private:
      //used to find new start node, if BFS cannot continue
      NodeId first_not_visited;
};
} // namespace partitioner