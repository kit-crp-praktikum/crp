#include "algorithms/dijkstra.hpp"
#include "crp/crp.h"
#include "partitioner/geo-data.h"
#include <memory>

namespace crp
{
void CRPAlgorithm::prepare(Graph *graph, partitioner::GeoData *geo_data)
{
    this->g = graph;
    this->bidir_dijkstra = std::make_unique<BidirectionalDijstkra>(g->num_nodes());
    this->partition = params.partitioner(graph, geo_data, params.number_of_levels + params.number_of_phantom_levels, params.cells_per_level);
}
} // namespace crp
