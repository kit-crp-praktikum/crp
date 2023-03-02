#include "cmdline-parsing.hpp"
#include "crp/crp.h"
#include "data-types.h"
#include "graph.h"
#include <fstream>
#include <iostream>
#include <map>
#include <queue>
#include <set>
#define _ << " " <<
#define debug(x) #x << " = " << x

struct PunchPartition
{
    int nr_cells;
    std::vector<crp::CellId> cell_for_node;

    void load_from_file(crp::Graph &g, std::string path, int level)
    {
        std::ifstream in(path);

        cell_for_node.assign(g.num_nodes(), -1);

        std::string unused;
        in >> unused >> nr_cells;

        std::vector<std::vector<NodeId>> border_nodes_of_cell(nr_cells);

        NodeId vertex;
        NodeId cellId;

        while (in >> unused >> vertex >> cellId)
        {
            border_nodes_of_cell[cellId].push_back(vertex);
            cell_for_node[vertex] = cellId;
        }

        for (auto cell = 0; cell < nr_cells; cell++)
        {
            for (auto x : border_nodes_of_cell[cell])
            {
                assign_cell_bfs(g, x, level);
            }
        }

        for (NodeId i = 0; i < g.num_nodes(); i++)
        {
            if (cell_for_node[i] == -1)
            {
                cell_for_node[i] = 0;
            }
        }
    }

  private:
    void assign_cell_bfs(crp::Graph &g, NodeId start, int level)
    {
        assert(cell_for_node[start] >= 0);
        std::queue<NodeId> q;
        q.push(start);

        int marked = 0;
        while (!q.empty())
        {
            auto x = q.front();
            q.pop();
            ++marked;

            for (auto [y, w] : g[x])
            {
                if (cell_for_node[y] == -1)
                {
                    cell_for_node[y] = cell_for_node[start];
                    q.push(y);
                }
            }
        }
    }
};

int main(int argc, char **argv)
{
    int pos;

    pos = find_required_argument(argc, argv, 'g', "graph", true);
    const std::string graph_path = argv[pos + 1];
    check_input_directory_valid(graph_path, "travel_time", false);

    pos = find_required_argument(argc, argv, 'p', "partition", true);
    const std::string partition_base = argv[pos + 1];

    pos = find_required_argument(argc, argv, 'l', "levels", true);
    const int levels = parse_integer_or_bail(argv[pos + 1]);

    auto g = load_graph_from_directory(graph_path, "travel_time");

    crp::RecursivePartitionMask mask(g.num_nodes());
    std::vector<PunchPartition> partitions(levels);

    // Load levels from the partition file
    for (int level = 1; level <= levels; level++)
    {
        partitions[level - 1].load_from_file(g, partition_base + "-level-" + std::to_string(level) + ".cut", level);
    }

    // Cell on a particular Level
    using Cell = std::pair<crp::CellId, crp::LevelId>;

    std::map<Cell, std::set<Cell>> children;

    for (int level = 0; level < levels; level++)
    {
        for (NodeId x = 0; x < g.num_nodes(); x++)
        {
            Cell cur = {partitions[level].cell_for_node[x], level};
            Cell upw = {(level + 1 < levels ? partitions[level + 1].cell_for_node[x] : 0), level + 1};
            children[upw].insert(cur);
        }
    }

    int cells_per_level = 1;
    for (auto &[cell, list] : children)
    {
        cells_per_level = std::max(cells_per_level, (int)list.size());
    }

    while (__builtin_popcount(cells_per_level) > 1)
        ++cells_per_level;
    int bits_per_level = 32 - __builtin_clz(cells_per_level - 1);
    std::cerr << "cells_per_level=" << cells_per_level << std::endl;
    std::cerr << "levels=" << levels << std::endl;
    std::cerr << "bits_per_level=" << bits_per_level << std::endl;

    std::map<Cell, int> nested_cell_id;
    for (auto &[x, list] : children)
    {
        int index = 0;
        for (auto &cell : list)
        {
            nested_cell_id[cell] = index;
            ++index;
        }
    }

    for (int level = 0; level < levels; level++)
    {
        for (NodeId i = 0; i < g.num_nodes(); i++)
        {
            mask[i] <<= bits_per_level;
            Cell cur = {partitions[level].cell_for_node[i], level};
            mask[i] |= nested_cell_id[cur];
        }
    }

    std::cout.write((char *)mask.data(), sizeof(uint32_t) * mask.size());
}
