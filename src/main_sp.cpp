
#include "algorithms/bellman_ford.hpp"
#include "algorithms/dijkstra.hpp"
#include "algorithms/dinics.hpp"
#include "algorithms/floyd_warshall.hpp"
#include "data-types.h"
#include "datastructure/timestamped_vector.hpp"
#include "graph.h"
#include "partitioner/inertial_flow.hpp"
#include "partitioner/preprocessing.hpp"
#include <cassert>
#include <iostream>
#include <limits>

void test_timestampled_v()
{
    int def = 1, n = 100;
    TimestampedVector<int> v(n, def);
    assert(v[0] == def);
    v[0] = 2;
    assert(v[0] == 2);
    v.reset();
    assert(v[0] == def);
}

void test_sp()
{
    int n = 10;
    using Node = std::pair<int, int>;
    std::vector<std::vector<Node>> gr(n);
    // path
    for (int i = 0; i < n - 1; i++)
    {
        gr[i].push_back({i + 1, 1});
        gr[i + 1].push_back({i, 1});
    }
    auto neighbors = [&](NodeId v, auto f) {
        for (auto [u, weight] : gr[v])
        {
            f(u, weight);
        }
    };
    int start = 5;
    Dijkstra dij(n);
    std::cout << "dijkstra \n";
    dij.compute_distance(start, neighbors);
    for (int i = 0; i < n; i++)
    {
        std::cout << dij.tentative_distance(i) << " ";
    }
    std::cout << "\n";

    BidirectionalDijstkra bdij(n);
    std::cout << "\nbidirectional dijkstra \n";
    for (int i = 0; i < n; i++)
    {
        auto [v, d] = bdij.compute_distance_target(start, i, neighbors, neighbors);
        std::cout << d << " ";
    }
    std::cout << "\n";

    BellmanFord bf(n);
    std::cout << "\nbellman ford \n";
    bf.compute_distance(start, neighbors);
    for (int i = 0; i < n; i++)
    {
        std::cout << bf.tentative_distance(i) << " ";
    }
    std::cout << "\n";

    FloydWarshall fw(n);
    std::cout << "\nfloyd warshall \n";
    fw.compute_all_distances(neighbors);
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            std::cout << fw.get_distance(i, j) << " ";
        }
        std::cout << "\n";
    }
}

void test_flow1()
{
    int n = 6;
    DinicsFlow dinics(n);
    dinics.add_edge(0, 1, 7);
    dinics.add_edge(0, 4, 4);
    dinics.add_edge(1, 2, 5);
    dinics.add_edge(1, 3, 3);
    dinics.add_edge(2, 5, 8);
    dinics.add_edge(3, 2, 3);
    dinics.add_edge(3, 5, 5);
    dinics.add_edge(4, 1, 3);
    dinics.add_edge(4, 3, 2);
    std::cout << "\ndinics"
              << "\n";
    std::cout << dinics.max_flow(0, n - 1) << "\n";
    auto [cut, v] = dinics.min_cut_partition(0, n - 1);
    for (auto b : v)
    {
        std::cout << b << " ";
    }
    std::cout << "\n";
}

void test_flow2()
{
    int n = 8;
    DinicsFlow dinics(n);
    // C_3 - 3 - 4 - C_3
    dinics.add_edge(0, 1, 1);
    dinics.add_edge(0, 2, 1);
    dinics.add_edge(1, 2, 1);
    dinics.add_edge(2, 3, 1);
    dinics.add_edge(3, 4, 1);
    dinics.add_edge(4, 5, 1);
    dinics.add_edge(5, 6, 1);
    dinics.add_edge(5, 7, 1);
    dinics.add_edge(6, 7, 1);
    std::vector<NodeId> src = {0, 1, 2};
    std::vector<NodeId> target = {5, 6, 7};
    std::cout << "\nmultistart dinics"
              << "\n";
    std::cout << dinics.multi_src_target_max_flow(src, target) << "\n";
    auto [cut, v] = dinics.multi_src_target_min_cut_partition(src, target);
    for (auto b : v)
    {
        std::cout << b << " ";
    }
    std::cout << "\n";
}

void test_inertial_flow()
{
    int n = 10;
    uint32_t lines = 10;
    double group_size = 0.1;
    using Node = std::pair<NodeId, Distance>;
    std::vector<std::vector<Node>> graph(n * n);
    std::vector<int> x(n * n);
    std::vector<int> y(n * n);
    auto add_edge = [&](NodeId v, NodeId w) {
        if (v < w)
        { // only in one direction
            graph[v].push_back({w, 1});
            graph[w].push_back({v, 1});
        }
    };
    int dx[4] = {0, 0, -1, 1};
    int dy[4] = {-1, 1, 0, 0};
    // grid graph
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            x[i * n + j] = i;
            y[i * n + j] = j;
            for (int k = 0; k < 4; k++)
            {
                int i1 = i + dx[k];
                int j1 = j + dy[k];
                if (std::min(i1, j1) >= 0 && std::max(i1, j1) < n)
                {
                    add_edge(i * n + j, i1 * n + j1);
                }
            }
        }
    }
    partitioner::GeoData geo_data(x, y);
    partitioner::InertialFlowPartitioner inertial(n * n, lines, group_size);
    std::cout << "\ninertial flow \n";
    std::vector<bool> v = inertial.partition(graph, geo_data);
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            std::cout << v[i * n + j];
        }
        std::cout << "\n";
    }
}

void test_preprocessing()
{
    /*
        0 -> 1    5 -- 6
        ^    |    |    |
        |    v    |    |
        3 <- 2 -> 4 -- 7
    */
    NodeId n = 10;
    using AdjacencyList = std::vector<std::vector<std::pair<NodeId, Distance>>>;
    AdjacencyList adj_list(n);
    Distance d = 1; // increasing edge weights
    std::vector<std::pair<NodeId, NodeId>> edges = {{0, 1}, {1, 2}, {2, 3}, {3, 0}, {2, 4}, {4, 5}, {5, 4},
                                                    {5, 6}, {6, 5}, {6, 7}, {7, 6}, {4, 7}, {7, 4}};
    for (auto [v, w] : edges)
    {
        adj_list[v].push_back({w, d++});
    }
    crp::Graph graph(adj_list);
    AdjacencyList new_graph = partitioner::make_undirected(graph);
    assert(new_graph.size() == adj_list.size());
    for (NodeId v = 0; v < n; v++)
    {
        for (auto [w, dist] : new_graph[v])
        {
            if (std::min(v, w) >= 4)
            {
                assert(dist == 2);
            }
            else
            {
                assert(dist == 1);
            }
            bool found = false;
            for (auto [u, dist2] : new_graph[w])
            {
                if (u == v && dist == dist2)
                {
                    found = true;
                    break;
                }
            }
            assert(found);
        }
    }
    std::cout << "\npreprocessing test ok \n";
}

int main()
{
    test_timestampled_v();
    test_sp();
    test_flow1();
    test_flow2();
    test_inertial_flow();
    test_preprocessing();
}