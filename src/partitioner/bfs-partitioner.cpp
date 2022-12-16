#pragma once

#include <vector>
#include <queue>
#include <algorithm>
#include "src/data-types.h"
#include "bfs-partitioner.h"
#include "src/graph.h"

partitioner::BfsPartitioner::BfsPartitioner(){}

std::vector<bool> partitioner::BfsPartitioner::partition(crp::Graph* g, std::vector<bool>* cont_nodes)
{


    int cnt_visited =1;
    int num_nodes = std::count((*cont_nodes).begin(), (*cont_nodes).end(), true);

    std::queue<NodeId> q;
    std::vector<bool> visited(g->num_nodes(), false);
    
    NodeId start = findStart(&visited, cont_nodes);
    visited[start] = true;
    q.push(start);

    while (cnt_visited < num_nodes/2)
    {
        NodeId s = q.front();
        q.pop();
        
        for (auto [to, w] : (*g)[s]) 
        { // edge s -> to with wieght w
            if (! (*cont_nodes)[to]) continue;
            if (! visited[to]) 
            {
                visited[to] = true;
                q.push(to);
                cnt_visited++;
            }
        }

        if (q.empty()) {
            // may occur if graph not connected
            start = findStart(&visited, cont_nodes);
            q.push(start);
            visited[start] = true;
        }
    }

    return visited;
}

NodeId partitioner::BfsPartitioner::findStart(std::vector<bool>* visited, std::vector<bool>* cont_nodes)
{
    int ind = 0;
    if (gd.latitude.empty()) {
        // find first unvisited
        while (!(*cont_nodes)[ind] || (*visited)[ind])
        {
            ind++;
        }
        return ind;
    }
    else 
    {
        // find highest latitude
        while (!(*cont_nodes)[ind]) ind++;

        int max_lat = gd.latitude[ind];
        for (int i = ind; i < gd.latitude.size(); i++) {
            if ((*cont_nodes)[ind] || gd.latitude[i] > max_lat) 
            {
                max_lat = gd.latitude[i];
                ind = i;
            }
        }
        return ind;
    }
}

