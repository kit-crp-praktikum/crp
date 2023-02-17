#pragma once

#include "algorithms/dinics.hpp"
#include "data-types.h"
#include "graph.h"
#include "partitioner/geo-data.h"
#include <algorithm>
#include <cassert>
#include <functional>
#include <math.h>
#include <numeric>
#include <vector>

#include "bipartitioner.h"
#include "lib/debug.h"

/**
 * Implementation of Inertial-Flow partitioner.
 */
namespace partitioner
{

struct InertialFlowParameters
{
    uint32_t number_of_lines;
    double group_size;
};
class InertialFlowPartitioner : public BiPartitioner
{
  public:
    /**
     * @param max_num_nodes maximum number of nodes which will be used during partitioning
     * @param number_of_lines number of lines which are used to sort the nodes. use each centered in the graph with
     * angles 2*PI/i i in {0, ..., num_lines - 1}
     * @param group_size size of the src/target node for the flow problem in percentage
     */
    InertialFlowPartitioner(uint32_t max_num_nodes, uint32_t number_of_lines, double group_size)
        : sorted_nodes(max_num_nodes), scalar_products(max_num_nodes), number_of_lines(number_of_lines),
          group_size(group_size)
    {
    }

    std::vector<bool> partition(crp::AdjacencyList &graph, partitioner::GeoData &geo_data)
    {
        uint32_t n = graph.size();
        assert(n == geo_data.latitude.size());
        assert(n == geo_data.longitude.size());
        assert(n <= sorted_nodes.size());
        assert(group_size < 0.5);
        assert(number_of_lines > 0);

        std::iota(sorted_nodes.begin(), sorted_nodes.begin() + n, 0);
        compute_center(geo_data);

        return minimize_edge_cut(graph, geo_data);
    }

    void set_number_of_lines(uint32_t num_lines)
    {
        number_of_lines = num_lines;
    }

    void set_group_size(double _group_size)
    {
        group_size = _group_size;
    }

  private:
    std::vector<bool> minimize_edge_cut(crp::AdjacencyList &graph, partitioner::GeoData &geo_data)
    {
        std::vector<NodeId> left_nodes, right_nodes;
        std::vector<bool> best_partition, partition;
        uint32_t best_balance = 0, balance;
        uint32_t best_cut = INF;

        uint32_t n = graph.size();
        DinicsFlow dinics(n);

        for (uint32_t i = 0; i < number_of_lines; i++)
        {
            double angle = i * (M_PI / number_of_lines);
            sort_by_line(cos(angle), sin(angle), geo_data);
            for (uint32_t i = 0; i < n && ((double)i / n) < group_size; i++)
            {
                left_nodes.push_back(sorted_nodes[i]);
                right_nodes.push_back(sorted_nodes[n - 1 - i]);
            }
            auto [cut_size, partition] = dinics.multi_src_target_min_cut_partition(graph, left_nodes, right_nodes);
            balance = std::count(partition.begin(), partition.end(), true);
            balance = std::min(balance, n - balance);
            if (cut_size < best_cut || ((cut_size == best_cut) && balance >= best_balance))
            {
                std::swap(best_partition, partition);
                std::swap(best_balance, balance);
                std::swap(best_cut, cut_size);
            }
            left_nodes.clear();
            right_nodes.clear();
        }
        return best_partition;
    }

    void compute_center(partitioner::GeoData &geo_data)
    {
        double n = geo_data.latitude.size();
        double center_x = std::reduce(geo_data.latitude.begin(), geo_data.latitude.end(), 0.0f, std::plus<float>());
        double center_y = std::reduce(geo_data.longitude.begin(), geo_data.longitude.end(), 0.0f, std::plus<float>());
        center_of_nodes = {center_x / n, center_y / n};
    }

    // line is defined by a * (x - center_x) + b * (y - center_y) = 0
    void sort_by_line(double a, double b, partitioner::GeoData &geo_data)
    {
        for (uint32_t i = 0; i < geo_data.latitude.size(); i++)
        {
            scalar_products[i] = a * (geo_data.latitude[i] - center_of_nodes.first) +
                                 b * (geo_data.longitude[i] - center_of_nodes.second);
        }
        auto compare = [&](const NodeId v, const NodeId w) { return scalar_products[v] < scalar_products[w]; };
        std::sort(sorted_nodes.begin(), sorted_nodes.begin() + geo_data.latitude.size(), compare);
    }

    std::vector<NodeId> sorted_nodes;
    std::vector<double> scalar_products;
    std::pair<double, double> center_of_nodes;
    uint32_t number_of_lines;
    double group_size;
};
} // namespace partitioner
