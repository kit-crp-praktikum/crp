#pragma once

#include "bfs-partitioner.h"

partitioner::BfsPartitioner::BfsPartitioner()
{
}

std::vector<bool> partitioner::BfsPartitioner::partition(crp::AdjacencyList &graph, GeoData &geo_data)
{
    // init
    uint32_t num_nodes = graph.size();
    uint32_t cnt_visited = 0;
    NodeId partition_size = num_nodes / 2;
    first_not_visited = 0;
    std::vector<bool> visited(num_nodes, false);
    std::queue<NodeId> q;

    NodeId start = findStart(visited, geo_data);
    visited[start] = true;
    q.push(start);
    cnt_visited++;

    while (cnt_visited < partition_size)
    {
        NodeId s = q.front();
        q.pop();

        for (auto [to, w] : graph[s])
        {
            // edge s -> to with weight w
            if (!visited[to])
            {
                visited[to] = true;
                q.push(to);
                cnt_visited++;
                if (cnt_visited >= partition_size)
                    break;
            }
        }

        if (cnt_visited < partition_size && q.empty())
        {
            // all component nodes visited, find another start node
            start = findStart(visited, geo_data);
            q.push(start);
            visited[start] = true;
            cnt_visited++;
        }
    }
    return visited;
}

NodeId partitioner::BfsPartitioner::findStart(std::vector<bool> &visited, GeoData &geo_data)
{
    // go to next unvisited node
    while (visited[first_not_visited])
        first_not_visited++;

    // find highest not visitied latitude node
    NodeId next_start = first_not_visited;
    float max_lat = geo_data.latitude[first_not_visited];
    for (unsigned i = first_not_visited + 1; i < geo_data.latitude.size(); i++)
    {
        if (!visited[i] && geo_data.latitude[i] > max_lat)
        {
            max_lat = geo_data.latitude[i];
            next_start = i;
        }
    }
    return next_start;
}
