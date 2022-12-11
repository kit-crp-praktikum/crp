
#include <iostream>
#include <cassert>
#include <limits>
#include "datastructure/timestamped_vector.hpp"
#include "algorithms/dijkstra.hpp"
#include "algorithms/bellman_ford.hpp"
#include "algorithms/floyd_warshall.hpp"
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


int main() {
    std::cout << "hi \n";

    test_timestampled_v();
    test_sp();

}