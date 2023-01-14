#pragma once

#include "geo-data.h"
#include "src/data-types.h"
#include "graph.h"
#include "geo-data.h"

#include <vector>

/**
 * Interface for BiPartitioner
 */
namespace partitioner
{
class BiPartitioner
{
  public:
    /**
     * Partitions the graph into two parts.
     */
    virtual std::vector<bool> partition(crp::AdjacencyList &graph, GeoData &geo_data) = 0;

    virtual ~BiPartitioner() = default;
};
} // namespace partitioner