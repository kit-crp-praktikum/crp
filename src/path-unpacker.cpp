#include "path-unpacker.h"
#include "data-types.h"
#include "graph.h"

crp::PathUnpackingResult crp::isPathCorrect(std::vector<NodeId> *path, Graph *g, Distance dist)
{
    NodeId u, v;
    Distance path_length = 0;
    for (unsigned i = 0; i < path->size() - 1; i++)
    {
        u = (*path)[i];
        v = (*path)[i + 1];

        auto w = g->get_edge(u, v);
        if (!w)
        {
            return PathUnpackingResult::EdgeMissing;
        }
        else
        {
            path_length += w.value();
        }
    }

    if (dist != path_length)
        return PathUnpackingResult::TotalLengthWrong;
    return PathUnpackingResult::Ok;
}
