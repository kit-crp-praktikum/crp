#pragma once

#include "algorithms/dijkstra.hpp"
#include "data-types.h"
#include "graph.h"
#include "partitioner/geo-data.h"
#include "shortest-path-algorithm.h"
#include <cstdint>
#include <functional>
#include <memory>
#include <span>

namespace crp
{

// A vector of bitmasks indicating the partitions each node is on.
using RecursivePartitionMask = std::vector<uint32_t>;

struct RecursivePartition
{
    int number_of_levels;
    int cells_per_level;
    RecursivePartitionMask mask;

    // Find the first level where u and v belong to different cells
    inline int find_level_differing(NodeId u, NodeId v)
    {
        const uint32_t diff = mask[u] ^ mask[v];
        if (!diff)
        {
            // On level -1 every node is separate
            return -1;
        }

        const uint32_t first_diff = __builtin_ctz(diff) / number_of_levels;
        return number_of_levels - first_diff - 1;
    }

    inline int find_cell_for_node(NodeId u, int level)
    {
        const uint32_t bits_per_level = get_bits_per_level();
        uint32_t bits = (number_of_levels - level) * bits_per_level;
        return mask[u] & ((1 << bits) - 1);
    }

    inline uint32_t get_bits_per_level() const
    {
        return 32 - __builtin_clz(cells_per_level - 1);
    }
};

// A function which can compute a partition from (graph, number_of_levels, cells_per_level).
using RecursivePartitionerFunction =
    std::function<RecursivePartition(crp::Graph *, partitioner::GeoData *data, int, int)>;

using LevelId = int;
using CellId = int;

struct OverlayStructure
{
    // Create an overlay structure using the given graph and a list of bitmasks for each node, indicating
    // which cell each node is placed on each level.
    OverlayStructure(crp::Graph *g, RecursivePartition partition);

    CellId get_cell_for_node(NodeId u, LevelId level);

    // Get the internal ID of the given border node inside its cell on the given level.
    NodeId get_internal_id(NodeId u, LevelId level);

    // Get a reference to the memory where the distance between border nodes a and b is stored.
    // Note that a and b are internal IDs for the given cell.
    Distance *get_distance(LevelId level, CellId cell, NodeId a, NodeId b);

    // Get the number of cells on a given level.
    // The cells are numbered from 0 to num_cells_in_level-1.
    int num_cells_in_level(int level);

    // Get a list of the border nodes of the given cell on the given level.
    // The i-th node in the list returned is the node with internal ID i.
    std::span<NodeId> get_border_nodes_for_cell(LevelId level, CellId cell);

    // Get a list of the cellIds contained in the current cell at level - 1
    std::span<CellId> get_child_cellIds(LevelId level, CellId cell);

    // Get a list of nodes in level 0 cell.
    std::span<NodeId> get_nodes_level0(CellId cell);

    // Returns the number of levels in the recursive partition.
    int get_number_of_levels();

  private:
    using Clique = std::vector<std::vector<Distance>>;

    // Store the distance between each pair of border nodes in each cell of each level.
    // Usage: cliques[level][cell]
    std::vector<std::vector<Clique>> cliques;

    // A list of NodeIDs for each node on each level. Used to figure out the index of a node in the matrices
    // representing cliques.
    std::vector<std::vector<NodeId>> node_id_on_level;

    // border_nodes[level][cell] contains a list of the border nodes of the given cell (on the given level).
    //
    // Cells on each level are numbered from 0 to the number of cells on that level.
    // Levels are numbered 0..(number_of_levels-1), where level -1 contains isolated nodes only.
    std::vector<std::vector<std::vector<NodeId>>> border_nodes;

    // child_cell_ids[level][cell]
    // contains all cellIds of level - 1, which are contained in the current cell
    std::vector<std::vector<std::vector<CellId>>> child_cell_ids;

    // nodes_in_level_0[cell]
    // contains all nodes in level 0 cell
    std::vector<std::vector<NodeId>> nodes_in_level_0;

    RecursivePartition partition;
};

// A function which runs the customization on the given graph
using CustomizationFunction = std::function<void(crp::Graph *, OverlayStructure *)>;

struct CRPAlgorithmParams
{
    int number_of_levels;
    int cells_per_level;
    RecursivePartitionerFunction partitioner;
    CustomizationFunction customizer;
};

/**
 * The actual implementation of the CRP Algorithm.
 */
class CRPAlgorithm : public CRPAlgorithmInterface
{
  public:
    void prepare(Graph *graph, partitioner::GeoData *geo_data) override;
    void customize() override;
    Distance query(NodeId start, NodeId end) override;
    std::vector<NodeId> query_path(NodeId start, NodeId end, Distance &out_dist) override;

  public:
    CRPAlgorithm(CRPAlgorithmParams params);

  private:
    CRPAlgorithmParams params;

    // Query the shortest path from start to end using bidir_dijkstra.
    // The result is the node where forward and backward search met and the computed shortest distance.
    template <bool update_parents> std::pair<NodeId, Distance> _query(NodeId start, NodeId end);

    // Finds path between start,end nodes on given level using bidir dijkstra
    // Inside their common cell on level+1
    std::vector<NodeId> _unpack (NodeId start, NodeId end, int level);

    crp::Graph *g;
    crp::Graph reverse;
    RecursivePartition partition;
    std::unique_ptr<OverlayStructure> overlay;
    std::unique_ptr<BidirectionalDijstkra> bidir_dijkstra;
};

// customization variants

// run dijkstra on original graph
void customize_with_dijkstra(crp::Graph *g, crp::OverlayStructure *overlay);

// run bellman-ford on original graph
void customize_with_bellman_ford(crp::Graph *g, crp::OverlayStructure *overlay);

// run dijkstra on rebuild cell-graph
void customize_dijkstra_rebuild(crp::Graph *g, crp::OverlayStructure *overlay);

// run bellman-ford on rebuild cell-graph
void customize_bellman_ford_rebuild(crp::Graph *g, crp::OverlayStructure *overlay);

// run floyd-warshall on rebuild cell-graph
void customize_floyd_warshall_rebuild(crp::Graph *g, crp::OverlayStructure *overlay);

} // namespace crp
