#pragma once

#include "data-types.h"
#include "graph.h"
#include "lib/timer.h"
#include <iostream>
#include <limits>
#include <queue>
#include <vector>

/**
 * Implementation of dinics-flow algorithm for undirected graphs.
 */
class DinicsFlow
{
    using EdgeId = uint32_t;
    using Flow = int32_t;

    struct Edge
    {
        NodeId to;
        EdgeId reverse_edge; // position in graph[to]
        Flow flow, capacity;
        Edge(NodeId _to, EdgeId _rev, Flow _cap) : to(_to), reverse_edge(_rev), flow(0), capacity(_cap)
        {
        }
        Flow remaining_capacity()
        {
            return capacity - flow;
        }
        bool has_capacity()
        {
            return remaining_capacity() > 0;
        }
    };

  public:
    DinicsFlow(uint32_t num_nodes)
        : graph(num_nodes + 2), // 2 extra nodes for multi src and target nodes
          work(num_nodes + 2), bfs_distance(num_nodes + 2), multi_src(num_nodes), multi_target(num_nodes + 1),
          is_source(num_nodes), is_target(num_nodes), num_nodes(num_nodes),
          INF_FLOW(std::numeric_limits<Flow>::max() / 2)
    {
    } // / 2 to avoid INF - flow overflow for flow < 0

    void add_edge(NodeId from, NodeId to, Flow capacity)
    {
        graph[from].emplace_back(to, graph[to].size(), capacity);
        graph[to].emplace_back(from, graph[from].size() - 1, capacity);
    }

    Flow max_flow(NodeId src, NodeId target)
    {
        reset();
        Flow total_flow = 0;
        Flow new_flow = 0;
        while (bfs(src, target))
        {
            std::fill(work.begin(), work.end(), 0);
            do
            {
                new_flow = dfs(src, target, INF_FLOW);
                total_flow += new_flow;
            } while (new_flow > 0);
        }
        return total_flow;
    }

    Flow multi_src_target_max_flow(const crp::AdjacencyList &original_graph, std::vector<NodeId> &src_nodes,
                                   std::vector<NodeId> &target_nodes)
    {
        setup_multi_src_target(original_graph, src_nodes, target_nodes);
        Flow flow = max_flow(multi_src, multi_target);
        reset_multi_src_target();
        return flow;
    }

    std::pair<uint32_t, std::vector<bool>> min_cut_partition(NodeId src, NodeId target)
    {
        Flow cut_size = max_flow(src, target);
        auto is_residual = [](Edge &e) { return e.has_capacity(); };
        std::vector<bool> left_most(num_nodes);
        generic_bfs(src, is_residual);
        for (NodeId v = 0; v < num_nodes; v++)
        {
            left_most[v] = (is_source[v] || bfs_distance[v] != -1);
        }
        return {cut_size, left_most};
    }

    std::pair<uint32_t, std::vector<bool>> multi_src_target_min_cut_partition(const crp::AdjacencyList &original,
                                                                              std::vector<NodeId> &src_nodes,
                                                                              std::vector<NodeId> &target_nodes)
    {
        setup_multi_src_target(original, src_nodes, target_nodes);
        auto [cut_size, partition] = min_cut_partition(multi_src, multi_target);

        reset_multi_src_target();
        return {cut_size, partition};
    }

  private:
    void setup_multi_src_target(const crp::AdjacencyList &original, std::vector<NodeId> &src_nodes,
                                std::vector<NodeId> &target_nodes)
    {
        for (NodeId v : src_nodes)
        {
            is_source[v] = 1;
        }

        for (NodeId v : target_nodes)
        {
            is_target[v] = 1;
        }

        for (NodeId u = 0; u < num_nodes; u++)
        {
            for (auto [v, w] : original[u])
            {
                if (is_source[u] && !is_source[v])
                {
                    add_edge(multi_src, v, 1);
                }

                if (!is_target[u] && is_target[v])
                {
                    add_edge(u, multi_target, 1);
                }

                if (!is_source[u] && !is_target[v] && u < v)
                {
                    add_edge(u, v, 1);
                }
            }
        }
    }

    void reset_multi_src_target()
    {
        std::fill(is_source.begin(), is_source.end(), false);
        std::fill(is_target.begin(), is_target.end(), false);
        graph.assign(num_nodes + 2, {});
    }

    void reset()
    {
        for (auto &v : graph)
        {
            for (Edge &e : v)
            {
                e.flow = 0;
            }
        }
    }

    void generic_bfs(NodeId start, auto predicate)
    {
        std::fill(bfs_distance.begin(), bfs_distance.end(), -1);
        queue.push(start);
        bfs_distance[start] = 0;
        while (!queue.empty())
        {
            NodeId v = queue.front();
            queue.pop();
            if (v == multi_target)
                continue;

            for (Edge &e : graph[v])
            {
                if (bfs_distance[e.to] == -1 && predicate(e))
                {
                    bfs_distance[e.to] = 1 + bfs_distance[v];
                    queue.push(e.to);
                }
            }
        }
    }

    // returns wether target is reachable
    bool bfs(NodeId src, NodeId target)
    {
        generic_bfs(src, [](Edge &e) { return e.has_capacity(); });
        return bfs_distance[target] != -1;
    }

    void push_flow(Edge &e, Flow flow)
    {
        e.flow += flow;
        graph[e.to][e.reverse_edge].flow -= flow;
    }

    // returns wether pushed flow
    Flow dfs(NodeId v, NodeId target, Flow residual_flow)
    {
        if (v == target || residual_flow == 0)
            return residual_flow;
        for (; work[v] < graph[v].size(); work[v]++)
        {
            Edge &e = graph[v][work[v]];
            if ((e.to == target || bfs_distance[e.to] == 1 + bfs_distance[v]) && e.has_capacity())
            {
                Flow flow = dfs(e.to, target, std::min(residual_flow, e.remaining_capacity()));
                if (flow > 0)
                {
                    push_flow(e, flow);
                    return flow;
                }
            }
        }
        return 0;
    }
    std::vector<std::vector<Edge>> graph;
    std::vector<NodeId> work; // for every node, what is the next neighbor to process
    std::vector<int32_t> bfs_distance;
    std::queue<NodeId> queue;
    NodeId multi_src;
    NodeId multi_target;
    std::vector<bool> is_source;
    std::vector<bool> is_target;
    NodeId num_nodes;
    Flow INF_FLOW;
};
