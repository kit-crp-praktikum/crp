#include "algorithms/dijkstra.hpp"
#include "crp/crp.h"
#include <memory>

namespace crp
{
void CRPAlgorithm::prepare(Graph *graph)
{
    this->g = graph;
    this->bidir_dijkstra = std::make_unique<BidirectionalDijstkra>(g->num_nodes());
    this->partition = params.partitioner(graph, params.number_of_levels, params.cells_per_level);
    this->overlay = std::make_unique<OverlayStructure>(graph, this->partition);
}
} // namespace crp
