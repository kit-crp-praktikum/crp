#pragma once

#include "data-types.h"
#include "graph.h"
#include <vector>

namespace crp
{
enum class PathUnpackingResult
{
    Ok = 0,
    EdgeMissing = 1,
    TotalLengthWrong = 2,
};

/**
 * Checks if found path is present in graph and length equals calculated distance
 */
PathUnpackingResult isPathCorrect(std::vector<NodeId> *path, Graph *g, Distance dist);
} // namespace crp
