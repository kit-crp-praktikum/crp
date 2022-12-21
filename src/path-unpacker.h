#pragma once

#include <vector>
#include "data-types.h"
#include "graph.h"

namespace crp
{
class PathUnpacker
{
  public:

    /**
     * Checks if found path is present in graph and length equals calculated distance
    */
    bool isPathCorrect (std::vector<NodeId>* path, Graph* g, Distance dist);

};
}