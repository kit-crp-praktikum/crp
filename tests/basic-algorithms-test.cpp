#include "partitioner/geo-data.h"
#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <doctest/doctest.h>

#include "algorithm-test.hpp"
#include "algorithms/bellman_ford.hpp"
#include "algorithms/bellman_ford_simd.hpp"
#include "algorithms/dijkstra.hpp"
#include "algorithms/floyd_warshall.hpp"

class BidirDijkstraAlgo : public crp::CRPAlgorithmInterface
{
  public:
    void prepare(crp::Graph *graph, partitioner::GeoData *)
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
        auto r = bi_dijkstra->compute_distance_target(
            start, end,
            [&](NodeId u, auto ForEach) {
                for (auto [v, weight] : (*g)[u])
                {
                    ForEach(v, weight);
                }
            },
            [&](NodeId u, auto ForEach) {
                for (auto [v, weight] : reversed[u])
                {
                    ForEach(v, weight);
                }
            });

        return r.second;
    }

    std::vector<NodeId> query_path(NodeId start, NodeId end, Distance &out_dist)
    {
        auto [middle, distance] = bi_dijkstra->compute_distance_target<true>(
            start, end,
            [&](NodeId u, auto ForEach) {
                for (auto [v, weight] : (*g)[u])
                {
                    ForEach(v, weight);
                }
            },
            [&](NodeId u, auto ForEach) {
                for (auto [v, weight] : reversed[u])
                {
                    ForEach(v, weight);
                }
            });
        return bi_dijkstra->unpack(start,end,middle);
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
    void prepare(crp::Graph *graph, partitioner::GeoData *)
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
            bf->compute_distance(s, [=](NodeId u, auto F) {
                for (auto [v, weight] : (*g)[u])
                {
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

    std::vector<NodeId> query_path(NodeId start, NodeId end, Distance &out_dist)
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
    void prepare(crp::Graph *graph, partitioner::GeoData *)
    {
        this->g = graph;
        this->fs = std::make_unique<FloydWarshall>(g->num_nodes());
    }

    void customize()
    {
        fs->compute_all_distances([=](NodeId u, auto F) {
            for (auto [v, weight] : (*g)[u])
            {
                F(v, weight);
            }
        });
    }

    Distance query(NodeId start, NodeId end)
    {
        return fs->get_distance(start, end);
    }

    std::vector<NodeId> query_path(NodeId start, NodeId end, Distance &out_dist)
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


class BellmanFordSIMDAlgo : public crp::CRPAlgorithmInterface
{
  public:
    void prepare(crp::Graph *graph, partitioner::GeoData *)
    {
        this->g = graph;
        this->bf = std::make_unique<BellmanFordSIMD>(g->num_nodes());
    }

    void customize()
    {
        int n = g->num_nodes();
        distance_matrix.resize(n, std::vector<Distance>(n));

        auto nghr = [&](NodeId u, auto F) {
            for (auto [v, weight] : (*g)[u])
            {
                F(v, weight);
            }
        };

        auto clamp = [&](int v) {
            return std::min(v, g->num_nodes() - 1);
        };
        std::array<NodeId, SIMD_LEN> start_nodes{};
        for (int s = 0; s < g->num_nodes(); s += SIMD_LEN)
        {
            for(int i = 0; i < SIMD_LEN; i++)
            {
                start_nodes[i] = clamp(s + i);
            }
            bf->compute_distance(start_nodes, nghr);
            for (int t = 0; t < g->num_nodes(); t++)
            {
                std::array<NodeId, SIMD_LEN> result = bf->tentative_distance(t);
                for(int i = 0; i < SIMD_LEN; i++)
                {   
                    distance_matrix[start_nodes[i]][t] = result[i];
                }
            }
        }
    }

    Distance query(NodeId start, NodeId end)
    {
        return distance_matrix[start][end];
    }

    std::vector<NodeId> query_path(NodeId start, NodeId end, Distance &out_dist)
    {
        // Not implemented yet
        return {};
    }

  private:
    std::unique_ptr<BellmanFordSIMD> bf;
    crp::Graph *g;
    std::vector<std::vector<Distance>> distance_matrix;
};

TEST_CASE("Test Bellman-Ford-SIMD")
{
    test_algorithm(std::make_unique<BellmanFordSIMDAlgo>());
}