#include "partitioner/geo-data.h"
#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <doctest/doctest.h>

#include "algorithm-test.hpp"
#include "crp/crp.h"
#include "grid-graph.hpp"
#include "partitioner/rec-partitioner.h"
#include <memory>

#include <iostream>
#define _ << " " <<
#define debug(x) #x << " = " << x

static crp::CRPAlgorithmParams gen_params_8x8()
{
    const int n = 8;
    crp::CRPAlgorithmParams params;
    params.number_of_levels = 2;
    params.cells_per_level = 4;
    params.partitioner = [](crp::Graph *g, partitioner::GeoData *geo, int nr_levels, int cells_per_level) {
        crp::RecursivePartition rp;
        rp.cells_per_level = cells_per_level;
        rp.number_of_levels = nr_levels;
        rp.mask = generate_two_level_partition(n);
        return rp;
    };

    params.customizer = crp::customize_with_dijkstra;

    dump_grid_graph_node_ids(n);
    return params;
}

TEST_CASE("Test CRP on 8x8 grid graph")
{
    auto crp = std::make_unique<crp::CRPAlgorithm>(gen_params_8x8());
    auto graph = generate_grid_graph(8);
    test_algorithm(std::move(crp), crp::Graph{graph});
}

TEST_CASE("Test CRP path unpacking on 8x8 grid graph")
{
    auto crp = std::make_unique<crp::CRPAlgorithm>(gen_params_8x8());
    auto graph = generate_grid_graph(8);
    test_algorithm_path(std::move(crp), crp::Graph{graph});
}
