#include "graph.h"
#include "data-types.h"
#include "lib/vector_io.h"

crp::Graph::Graph(std::string first_out, std::string head, std::string weights)
{
    this->first_out = load_vector<uint32_t>(first_out);
    this->head = load_vector<NodeId>(head);
    this->weights = load_vector<Distance>(weights);
}
crp::Graph::Graph(
    const std::vector<std::vector<std::pair<NodeId, Distance>>>& adj_list)
{
    uint32_t n = adj_list.size();
    for (uint32_t i = 0; i < n; i++)
    {
        first_out.push_back(head.size());
        for (auto [to, w] : adj_list[i]) {
            head.push_back(to);
            weights.push_back(w);
        }
    }

    // Sentinel at the end
    first_out.push_back(head.size());
}
