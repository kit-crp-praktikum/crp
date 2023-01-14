#include "data-types.h"
#include "graph.h"
#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <doctest/doctest.h>
#include <iomanip>
#include <numeric>

#include "crp/crp.h"
#include "grid-graph.hpp"

TEST_CASE("Grid graph n=8")
{
    const int n = 8;
    crp::Graph g{generate_grid_graph(n)};

    dump_grid_graph_node_ids(n);

    // Generate simple 2-level partition
    crp::RecursivePartition rp;
    rp.cells_per_level = 4;
    rp.number_of_levels = 2;
    rp.mask = generate_two_level_partition(n);

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
    // As the level 0 cells are 2x2, all of the nodes are border nodes to other cells, except
    // the nodes in the corners.
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
                    if (cur != 0 && cur != 7 && cur != 56 && cur != 63)
                    {
                        border_nodes0.back().push_back(cur);
                    }
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
