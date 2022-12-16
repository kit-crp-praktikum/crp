#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <doctest/doctest.h>

#include "algorithms/dijkstra.hpp"
#include "data-types.h"
#include <memory>
#include "graph.h"
#include "shortest-path-algorithm.h"
#include "graph-generator.hpp"

void test_algorithm(std::unique_ptr<crp::CRPAlgorithm> algorithm)
{
    crp::Graph g = generate_random_undirected_graph(100, 300, 100);

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


    algorithm->prepare(&g);
    algorithm->customize();

    for (int a = 0; a < g.num_nodes(); a++) {
        for (int b = 0; b < g.num_nodes(); b++) {
            CHECK(algorithm->query(a, b) == shortest_path_dijkstra(a, b));
        }
    }
}

class BidirDijkstraAlgo : public crp::CRPAlgorithm
{
  public:
    void prepare(crp::Graph *graph)
    {
        this->g = graph;
        this->reversed = graph->reversed();
        bi_dijkstra = std::make_unique<BidirectionalDijstkra>(graph->num_nodes());
    }

    void customize()
    {
        // Nothing to do here for dijkstra
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
