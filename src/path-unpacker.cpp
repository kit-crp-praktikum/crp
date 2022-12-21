#include "graph.h"
#include "data-types.h"
#include "path-unpacker.h"

bool crp::PathUnpacker::isPathCorrect(std::vector<NodeId>* path, Graph* g, Distance dist) {
    NodeId u,v;
    Distance path_length = 0;
    for (unsigned i = 0; i < path->size()-1; i++) {
        u = (*path)[i];
        v = (*path)[i+1];

        // check if {u,v} edge present in graph
        bool present = false;
        for (auto [to, w] : (*g)[u]) { // edge i -> to with wieght w
            if (to == v) {
                present = true; 
                path_length += w;
                break;
            }
        }
        if (!present) return false;
    }
    
    if (dist != path_length) return false;
    return true;
}