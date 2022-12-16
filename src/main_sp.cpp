
#include <iostream>
#include <cassert>
#include <limits>
#include "datastructure/timestamped_vector.hpp"
#include "algorithms/dijkstra.hpp"
#include "algorithms/bellman_ford.hpp"
#include "algorithms/floyd_warshall.hpp"
#include "algorithms/dinics.hpp"
#include "data-types.h"

void test_timestampled_v() {
    int def = 1, n = 100;
    TimestampedVector<int> v(n, def);
    assert(v[0] == def);
    v[0] = 2;
    assert(v[0] == 2);
    v.reset();
    assert(v[0] == def);
}

void test_sp() {
    int n = 10;
    using Node = std::pair<int, int>;
    std::vector<std::vector<Node>> gr(n);
    //path
    for(int i = 0; i < n - 1; i++) {
        gr[i].push_back({i + 1, 1});
        gr[i + 1].push_back({i, 1});
    }
    auto neighbors = [&](NodeId v, auto f) {
        for(auto [u, weight] : gr[v]) {
            f(u, weight);
        }
    };
    int start = 5;
    Dijkstra dij(n);
    std::cout << "dijkstra \n";
    dij.compute_distance(start, neighbors);
    for(int i = 0; i < n; i++) {
        std::cout << dij.tentative_distance(i) << " ";
    }
    std::cout << "\n";


    BidirectionalDijstkra bdij(n);
    std::cout << "\nbidirectional dijkstra \n";
    for(int i = 0; i < n; i++) {
        auto[v, d] = bdij.compute_distance_target(start, i, neighbors, neighbors);
        std::cout << d << " ";
    }
    std::cout << "\n";


    BellmanFord bf(n);
    std::cout << "\nbellman ford \n";
    bf.compute_distance(start, neighbors);
    for(int i = 0; i < n; i++) {
        std::cout << bf.tentative_distance(i) << " ";
    }
    std::cout << "\n";

    FloydWarshall fw(n);
    std::cout << "\nfloyd warshall \n";
    fw.compute_all_distances(neighbors);
    for(int i = 0; i < n; i++) {
        for(int j = 0; j < n; j++) {
            std::cout << fw.get_distance(i, j) << " ";
        }
        std::cout << "\n";
    }
}

void test_flow1() {
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
    std::cout << "\ndinics" << "\n";
    std::cout << dinics.max_flow(0, n - 1) << "\n";
    auto v = dinics.min_cut_partition(0, n - 1);
    for(auto b : v) {
        std::cout << b << " ";
    }
    std::cout << "\n";  
}

void test_flow2() {
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
    std::cout << "\nmultistart dinics" << "\n";
    std::cout << dinics.multi_src_target_max_flow(src, target) << "\n";
    auto v = dinics.multi_src_target_min_cut_partition(src, target);
    for(auto b : v) {
        std::cout << b << " ";
    }
    std::cout << "\n";  

}


int main() {
    std::cout << "hi \n";

    test_timestampled_v();
    test_sp();
    test_flow1();
    test_flow2();

}