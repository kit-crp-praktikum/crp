#pragma once

#include "bfs-partitioner.h"
#include "src/data-types.h"
#include <algorithm>
#include <queue>
#include <unordered_set>
#include <vector>

partitioner::BfsPartitioner::BfsPartitioner()
{
}

std::pair<std::unordered_set<NodeId> *, std::unordered_set<NodeId> *> partitioner::BfsPartitioner::partition(
    std::vector<std::vector<std::pair<NodeId, Distance>>> *g)
{
    int cnt_visited = 1;
    int num_nodes = g->size();
    std::queue<NodeId> q;
    std::vector<bool> visited(g->size(), false);

    result0.clear();
    result1.clear();

    NodeId start = findStart(&visited);
    visited[start] = true;
    q.push(start);

    while (cnt_visited < num_nodes / 2)
    {
        NodeId s = q.front();
        q.pop();

        for (auto [to, w] : (*g)[s])
        { // edge s -> to with wieght w
            if (!visited[to])
            {
                visited[to] = true;
                q.push(to);
                cnt_visited++;
                if (cnt_visited >= num_nodes / 2)
                    break;
            }
        }

        if (cnt_visited < num_nodes / 2 && q.empty())
        {
            // all component nodes visited, find another start node
            start = findStart(&visited);
            q.push(start);
            visited[start] = true;
            cnt_visited++;
        }
    }

    for (unsigned i = 0; i < visited.size(); i++)
        visited[i] ? result1.insert(i) : result0.insert(i);

    return {&result0, &result1};
}

NodeId partitioner::BfsPartitioner::findStart(std::vector<bool> *visited)
{
    if (gd.latitude.empty())
    {
        // find first unvisited
        while ((*visited)[lastStart])
        {
            lastStart++;
        }
        return lastStart;
    }
    else
    {
        // find highest latitude
        while ((*visited)[lastStart])
            lastStart++;

        int max_lat = gd.latitude[lastStart];
        for (unsigned i = lastStart + 1; i < gd.latitude.size(); i++)
        {
            if (!(*visited)[i] && gd.latitude[i] > max_lat)
            {
                max_lat = gd.latitude[i];
                lastStart = i;
            }
        }
        return lastStart;
    }
}
