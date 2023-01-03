#include "data-types.h"
#include "graph.h"
#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <doctest/doctest.h>
#include <iomanip>
#include <numeric>

#include "crp/crp.h"

#define encode(i, j) (n * (i) + j)

// Generate an n x n grid graph with all equal weights
crp::AdjacencyList generate_grid_graph(int n)
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

TEST_CASE("Grid graph n=8")
{
    const int n = 8;
    crp::Graph g{generate_grid_graph(n)};

    // Generate simple 2-level partition
    crp::RecursivePartition rp;
    rp.cells_per_level = 4;
    rp.number_of_levels = 2;
    rp.mask.resize(g.num_nodes());
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            std::cout << std::setw(2) << std::setfill(' ') << encode(i, j) << " ";
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
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
            rp.mask[cur] = mask;
            std::cout << std::setw(2) << std::setfill(' ') << mask << " ";
        }
        std::cout << std::endl;
    }

    crp::OverlayStructure os(&g, rp);

    CHECK(os.num_cells_in_level(0) == 16);
    CHECK(os.num_cells_in_level(1) == 4);

    std::vector<std::vector<NodeId>> border_nodes1 = {
        {3, 11, 19, 24, 25, 26, 27},
        {4, 12, 20, 28, 29, 30, 31},
        {32, 33, 34, 35, 43, 51, 59},
        {36, 37, 38, 39, 44, 52, 60},
    };

    // Generate the border nodes for the level 0.
    // As the level 0 cells are 2x2, most of the nodes are border nodes
    // at that level.
    std::vector<std::vector<NodeId>> border_nodes0;
    for (int cell = 0; cell < 16; cell++)
    {
        border_nodes0.push_back({});
        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < n; j++)
            {
                int cur = encode(i, j);
                if (rp.mask[cur] == cell)
                {
                    if ((i % 4 == 0 || i % 4 == 3) && (j % 4 == 0 || j % 4 == 3))
                    {
                        continue;
                    }

                    border_nodes0.back().push_back(cur);
                }
            }
        }
    }

    std::vector<std::vector<std::vector<NodeId>>> border_nodes = {
        border_nodes0,
        border_nodes1,
    };

    for (int level = 0; level < 2; level++)
    {
        for (int cellId = 0; cellId < os.num_cells_in_level(level); cellId++)
        {
            auto bnodes = os.get_border_nodes_for_cell(level, cellId);
            REQUIRE(bnodes.size() == border_nodes[level][cellId].size());
            for (size_t i = 0; i < bnodes.size(); i++)
            {
                CHECK(bnodes[i] == border_nodes[level][cellId][i]);
            }
        }
    }

    for (int level = 0; level < 2; level++)
    {
        for (int cellId = 0; cellId < os.num_cells_in_level(level); cellId++)
        {
            auto bnodes = os.get_border_nodes_for_cell(level, cellId);
            REQUIRE(bnodes.size() == border_nodes[level][cellId].size());
            for (size_t i = 0; i < bnodes.size(); i++)
            {
                CHECK(bnodes[i] == border_nodes[level][cellId][i]);
            }
        }
    }

    // Check internal IDs
    for (int level = 0; level < 2; level++)
    {
        for (int cellId = 0; cellId < os.num_cells_in_level(level); cellId++)
        {
            std::vector<NodeId> nodes_in_cell = border_nodes[level][cellId];

            // After getting the internal IDs, they should be a permutation of 0 .. (num_nodes_in_cell - 1)
            for (auto &x : nodes_in_cell)
            {
                CHECK(os.get_cell_for_node(x, level) == cellId);
                x = os.get_internal_id(x, level);
            }

            std::sort(nodes_in_cell.begin(), nodes_in_cell.end());

            std::vector<NodeId> expected(nodes_in_cell.size());
            std::iota(expected.begin(), expected.end(), 0);
            CHECK(nodes_in_cell == expected);
        }
    }
}
