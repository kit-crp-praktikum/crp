#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <doctest/doctest.h>

#include "algorithms/dijkstra.hpp"
#include "algorithms/bellman_ford.hpp"
#include "algorithms/floyd_warshall.hpp"
#include "data-types.h"
#include <memory>
#include "graph.h"
#include "shortest-path-algorithm.h"
#include "graph-generator.hpp"

void test_algorithm(std::unique_ptr<crp::CRPAlgorithmInterface> algorithm)
{
    // First, generate a graph with edges only equal to 1.
    // After that, change the weights for the customization phase
    crp::Graph g = generate_random_undirected_graph(100, 300, 1);
    algorithm->prepare(&g);

    std::mt19937 mt(0);
    // Random integer in [l, r]
    const auto& rnd_integer = [&] (int l, int r) {
        return std::uniform_int_distribution<int>(l,r)(mt);
    };
    for (auto& x : g.weights) {
        x = rnd_integer(1, 100);
    }

    algorithm->customize();

    Dijkstra plain_dijkstra(g.num_nodes());
    const auto& shortest_path_dijkstra = [&] (NodeId a, NodeId b) {
        plain_dijkstra.reset();
        plain_dijkstra.compute_distance_target(a, b, [&] (NodeId u, auto ForEachNeighbor) {
            for (auto [v, weight] : g[u]) {
                ForEachNeighbor(v, weight);
            }
        });

        return plain_dijkstra.tentative_distance(b);
    };

    for (int a = 0; a < g.num_nodes(); a++) {
        for (int b = 0; b < g.num_nodes(); b++) {
            CHECK(algorithm->query(a, b) == shortest_path_dijkstra(a, b));
        }
    }
}

class BidirDijkstraAlgo : public crp::CRPAlgorithmInterface
{
  public:
    void prepare(crp::Graph *graph)
    {
        this->g = graph;
        bi_dijkstra = std::make_unique<BidirectionalDijstkra>(graph->num_nodes());
    }

    void customize()
    {
        this->reversed = g->reversed();
    }

    Distance query(NodeId start, NodeId end)
    {
        auto r = bi_dijkstra->compute_distance_target(start, end,
            [&] (NodeId u, auto ForEach) {
                for (auto [v, weight] : (*g)[u]) {
                    ForEach(v, weight);
                }
            },
            [&] (NodeId u, auto ForEach) {
                for (auto [v, weight] : reversed[u]) {
                    ForEach(v, weight);
                }
            });

        return r.second;
    }

    std::vector<NodeId> query_path(NodeId start, NodeId end, Distance& out_dist)
    {
        // Not implemented yet
        return {};
    }

  private:
    std::unique_ptr<BidirectionalDijstkra> bi_dijkstra;
    crp::Graph *g;
    crp::Graph reversed;
};

TEST_CASE("Test Bidirectional Dijkstra")
{
    test_algorithm(std::make_unique<BidirDijkstraAlgo>());
}

class BellmanFordAlgo : public crp::CRPAlgorithmInterface
{
  public:
    void prepare(crp::Graph *graph)
    {
        this->g = graph;
        this->bf = std::make_unique<BellmanFord>(g->num_nodes());
    }

    void customize()
    {
        int n = g->num_nodes();
        distance_matrix.resize(n, std::vector<Distance>(n));

        for (int s = 0; s < g->num_nodes(); s++)
        {
            bf->compute_distance(s, [=] (NodeId u, auto F)
            {
                for (auto [v, weight] : (*g)[u]) {
                    F(v, weight);
                }
            });

            for (int t = 0; t < g->num_nodes(); t++)
            {
                distance_matrix[s][t] = bf->tentative_distance(t);
            }
        }
    }

    Distance query(NodeId start, NodeId end)
    {
        return distance_matrix[start][end];
    }

    std::vector<NodeId> query_path(NodeId start, NodeId end, Distance& out_dist)
    {
        // Not implemented yet
        return {};
    }

  private:
    std::unique_ptr<BellmanFord> bf;
    crp::Graph *g;
    std::vector<std::vector<Distance>> distance_matrix;
};

TEST_CASE("Test Bellman-Ford")
{
    test_algorithm(std::make_unique<BellmanFordAlgo>());
}

class FloydWarshallAlgo : public crp::CRPAlgorithmInterface
{
  public:
    void prepare(crp::Graph *graph)
    {
        this->g = graph;
        this->fs = std::make_unique<FloydWarshall>(g->num_nodes());
    }

    void customize()
    {
        fs->compute_all_distances([=] (NodeId u, auto F)
        {
            for (auto [v, weight] : (*g)[u]) {
                F(v, weight);
            }
        });
    }

    Distance query(NodeId start, NodeId end)
    {
        return fs->get_distance(start, end);
    }

    std::vector<NodeId> query_path(NodeId start, NodeId end, Distance& out_dist)
    {
        // Not implemented yet
        return {};
    }

  private:
    std::unique_ptr<FloydWarshall> fs;
    crp::Graph *g;
};

TEST_CASE("Test Floyd-Warshall")
{
    test_algorithm(std::make_unique<FloydWarshallAlgo>());
}
