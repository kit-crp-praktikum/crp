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

TEST_CASE("Test CRP on 8x8 grid graph")
{
    const int n = 8;

    crp::CRPAlgorithmParams params;
    params.number_of_levels = 2;
    params.cells_per_level = 4;
    params.partitioner = [](crp::Graph *g, int nr_levels, int cells_per_level) {
        crp::RecursivePartition rp;
        rp.cells_per_level = cells_per_level;
        rp.number_of_levels = nr_levels;
        rp.mask = generate_two_level_partition(n);
        return rp;
    };

    params.customizer = [](crp::Graph *g, crp::OverlayStructure *os) {
        // nothing as of yet, since customizer is hardcoded
    };

    dump_grid_graph_node_ids(n);

    auto crp = std::make_unique<crp::CRPAlgorithm>(params);
    auto graph = generate_grid_graph(n);
    test_algorithm(std::move(crp), crp::Graph{graph});
}
