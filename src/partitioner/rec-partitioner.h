#pragma once

#include <vector>
#include <unordered_set>
#include <map>
#include "src/data-types.h"
#include "bfs-partitioner.h"

/**
 * Recursive partitioning using 2-partitioner interface
 * Bit mask for each node is created 
 * Lowest lewel bits are being pushed to the end of the mask
*/

using ClusterId = uint32_t;

namespace partitioner
{
template<class BiPartitioner>
class RecPartitioner {
    public:
    RecPartitioner(BiPartitioner* bip) : biPartitioner(bip) {}

    /**
     * Partition graph recursively
     * Create bit masks for partitions
     * @param cells_per_level: # of cells on each level
     * @param levels: # of levels
     * @param g: graph to be partitioned
     * @return vector of bit masks for each node
    */
    std::vector<NodeId>* partition_rec(
        int cells_per_level, int levels, std::vector<std::vector<std::pair<NodeId, Distance>>>* g)
    {
        masks.clear();
        masks.resize(g->size());

        // shift level info into masks
        bits_for_mask = 0;
        int shifted = cells_per_level;
        while (shifted > 0) {
            shifted = (shifted >> 1);
            bits_for_mask++;
        }

        partition_rec_impl(cells_per_level, levels, g);
        return &masks;
    }


private:
    /**
     * Returns a bitmask for each node indicating which cells it belongs to.
     * Private version for recursive calls
     * @param cells_per_level: # of cells on each level
     * @param levels: # of levels
     * @param g: graph to be partitioned
    */
    void partition_rec_impl(
        int cells_per_level, int levels, std::vector<std::vector<std::pair<NodeId, Distance>>>* g)
    {
        if (levels == 0) return;

        // variables for graphs and mappings of the new partitiones
        std::map<NodeId,NodeId> map0[cells_per_level], map1[cells_per_level];
        std::vector<std::vector<std::pair<NodeId, Distance>>> cell0[cells_per_level], cell1[cells_per_level];

        // contains each current cell's graph on this level, and node id mapping to original id-s
        std::vector< std::vector<std::vector<std::pair<NodeId, Distance>>>* > cells;
        std::vector<std::map<NodeId,NodeId>*> mappings;

        std::map<NodeId,NodeId> map_default;
        for (NodeId i = 0; i < g->size(); i++) {
            map_default.insert({i,i});
        }
        mappings.push_back(&map_default);
        cells.push_back(g);

        // divide largest cell in two until there are cells_per_level cells
        for (int cell_cnt = 1; cell_cnt < cells_per_level; cell_cnt++) {
            // find largest cell
            int divide_pos = 0;
            for (unsigned i = 1; i < cells.size(); i++) {
                if (cells[i]->size() > cells[divide_pos]->size()) divide_pos = i;
            }
            std::vector<std::vector<std::pair<NodeId, Distance>>>* cell_to_divide = cells[divide_pos];

            // partition largest cell
            auto [nodes0, nodes1] = (*biPartitioner).partition(cell_to_divide);

            // calculate new cells and mappings
            divide_graph(mappings[divide_pos], nodes0, nodes1, cell_to_divide,
                &map0[cell_cnt], &map1[cell_cnt], &cell0[cell_cnt], &cell1[cell_cnt]);
            
            // remove old cell, add two new cells
            cells.erase(cells.begin()+divide_pos);
            mappings.erase(mappings.begin()+divide_pos);
            cells.push_back(&cell0[cell_cnt]);
            mappings.push_back(&map0[cell_cnt]);
            cells.push_back(&cell1[cell_cnt]);
            mappings.push_back(&map1[cell_cnt]);
        }

        // update masks for current layer's cells
        for (ClusterId part_id = 0; part_id < mappings.size(); part_id++) {
            std::map<NodeId,NodeId>* part_mapping = mappings[part_id];
            for (auto [v_part,v_orig] : (*part_mapping)) {
                // vv: original index in graph
                masks[v_orig] = masks[v_orig] << bits_for_mask;
                masks[v_orig] = masks[v_orig] | part_id;
            }
        }        

        // divide cells into cells on lower level recursively
        for (int i = 0; i < cells_per_level; i++) {
            partition_rec_impl(cells_per_level, levels-1, cells[i]);
        }
        return;
    }


    private:
    /**
     * divide current graph into two smaller ones
     * there are three kinds of graphs, each with its own node labeling: 
     * - original: the graph at the very beginning
     * - actual: graph which is being divided
     * - new: the newly created two graphs
     * @param map_act_orig: maps actual graph indices to original ones
     * @param nodes0: nodes of first partition
     * @param nodes1: nodes of second partition
     * @param g: graph to divide (actual graph)
     * @return multiple values by pointers: map_new_orig_0, map_new_orig_1, g0, g1
     * 
    */
    void divide_graph (
        std::map<NodeId,NodeId>* map_act_orig, std::unordered_set<NodeId>* nodes0, std::unordered_set<NodeId>* nodes1,
        std::vector<std::vector<std::pair<NodeId, Distance>>>* g,
        std::map<NodeId,NodeId>* map_new_orig_0, std::map<NodeId,NodeId>* map_new_orig_1,
        std::vector<std::vector<std::pair<NodeId, Distance>>>* g0, std::vector<std::vector<std::pair<NodeId, Distance>>>* g1)
    {
        // create new_orig mapping
        std::map<NodeId,NodeId> map_act_new_0, map_act_new_1;
        NodeId new_node_id = 0;
        for (auto node : (*nodes0)) {
            //map_new_act_0.push_back(node);
            map_new_orig_0->insert({new_node_id, (*map_act_orig)[node]});
            map_act_new_0.insert({node, new_node_id});
            new_node_id++;
        } 
        new_node_id = 0;
        for (auto node : (*nodes1)) {
            map_new_orig_1->insert({new_node_id, (*map_act_orig)[node]});
            map_act_new_1.insert({node, new_node_id});
            new_node_id++;
        } 

        // create new graphs
        // find u->v edge in actual graph, add it to new graph
        g0->resize(nodes0->size());
        for (auto u_act : (*nodes0)) {
            NodeId u_new = map_act_new_0[u_act];
            for (unsigned i = 0; i < (*g)[u_act].size(); i++) {
                auto [v_act,weight] = (*g)[u_act][i];
                NodeId v_new = map_act_new_0[v_act];
                // if u,v come to the same part after division, add edge to new graph
                if (map_act_new_0.contains(v_act)) {
                    (*g0)[u_new].push_back({v_new,weight});
                }
            }
        }
        g1->resize(nodes1->size());
        for (auto u_act : (*nodes1)) {
            NodeId u_new = map_act_new_1[u_act];
            for (unsigned i = 0; i < (*g)[u_act].size(); i++) {
                auto [v_act,weight] = (*g)[u_act][i];
                NodeId v_new = map_act_new_1[v_act];
                // if u,v come to the same part after division, add edge to new graph
                if (map_act_new_1.contains(v_act)) {
                    (*g1)[u_new].push_back({v_new,weight});
                }
            }
        }
    }

    private:
    BiPartitioner* biPartitioner;
    std::vector<ClusterId> masks;
    int bits_for_mask;
};
}