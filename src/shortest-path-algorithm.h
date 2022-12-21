#pragma once

#include "data-types.h"
#include "graph.h"
#include <memory>

namespace crp
{
class CRPAlgorithmInterface
{
  public:
    virtual ~CRPAlgorithmInterface() = default;

    /**
     * Phase 1: Initialize the routing algorithm for the given graph.
     * Note that the graph's weights do not have to have been initialized at this point yet,
     * and can all be equal zero.
     */
    virtual void prepare(Graph *graph) = 0;

    /**
     * Phase 2: Signal to the algorithm that the graph's weights have changed and it should update its
     * internal state to reflect these changes.
     */
    virtual void customize() = 0;

    /**
     * Query: Compute the shortest distance between start and end.
     */
    virtual Distance query(NodeId start, NodeId end) = 0;

    /**
     * Query: Find the shortest path between start and end.
     */
    virtual std::vector<NodeId> query_path(NodeId start, NodeId end, Distance &out_dist) = 0;
};

}; // namespace crp
