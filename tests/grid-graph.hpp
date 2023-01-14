#pragma once

#include "crp/crp.h"
#include "graph.h"
#include "partitioner/geo-data.h"
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

// Generate embedding for n x n grid graph
inline partitioner::GeoData generate_grid_graph_embedding(int n)
{
    uint32_t nr_cells = n * n;
    partitioner::GeoData geo_data;
    geo_data.latitude.resize(nr_cells);
    geo_data.longitude.resize(nr_cells);

    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            geo_data.latitude[encode(i, j)] = i;
            geo_data.longitude[encode(i, j)] = j;
        }
    }
    return geo_data;
}

// N must be a power of two.
inline crp::RecursivePartitionMask generate_two_level_partition(int n)
{
    crp::RecursivePartitionMask result(n * n);
    int nr_bits = __builtin_ctz(n);
    const int level0_bit = (1 << (nr_bits - 2));
    const int level1_bit = (1 << (nr_bits - 1));

    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            int cur = encode(i, j);

            int mask = 0;
            mask |= !!(i & level0_bit);
            mask <<= 1;
            mask |= !!(j & level0_bit);
            mask <<= 1;

            mask |= !!(i & level1_bit);
            mask <<= 1;
            mask |= !!(j & level1_bit);
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
