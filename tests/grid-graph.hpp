#pragma once

#include "crp/crp.h"
#include "graph.h"
#include <iomanip>
#include <iostream>

#define encode(i, j) (n * (i) + j)

// Generate an n x n grid graph with all equal weights
inline crp::AdjacencyList generate_grid_graph(int n)
{
    int nr_cells = n * n;
    crp::AdjacencyList adj_list(nr_cells);

    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            int cur = encode(i, j);
            if (i > 0)
                adj_list[cur].emplace_back(encode(i - 1, j), 1);
            if (j > 0)
                adj_list[cur].emplace_back(encode(i, j - 1), 1);
            if (i < n - 1)
                adj_list[cur].emplace_back(encode(i + 1, j), 1);
            if (j < n - 1)
                adj_list[cur].emplace_back(encode(i, j + 1), 1);
        }
    }

    return adj_list;
}

inline crp::RecursivePartitionMask generate_two_level_partition_for_8x8()
{
    const int n = 8;
    crp::RecursivePartitionMask result(n * n);

    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            int cur = encode(i, j);

            int mask = 0;
            mask |= !!(i & 2);
            mask <<= 1;
            mask |= !!(j & 2);
            mask <<= 1;

            mask |= !!(i & 4);
            mask <<= 1;
            mask |= !!(j & 4);
            result[cur] = mask;
            std::cout << std::setw(2) << std::setfill(' ') << mask << " ";
        }
        std::cout << std::endl;
    }

    return result;
}

inline void dump_grid_graph_node_ids(int n)
{
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            std::cout << std::setw(2) << std::setfill(' ') << encode(i, j) << " ";
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
}
