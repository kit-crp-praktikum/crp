#include "graph.h"
#include "data-types.h"
#include "lib/vector_io.h"

crp::Graph::Graph(std::string first_out, std::string head, std::string weights)
{
    this->first_out = load_vector<uint32_t>(first_out);
    this->head = load_vector<NodeId>(head);
    this->weights = load_vector<Distance>(weights);
}
crp::Graph::Graph(const AdjacencyList& adj_list)
{
    NodeId n = adj_list.size();
    for (NodeId i = 0; i < n; i++) {
        first_out.push_back(head.size());
        for (auto [to, w] : adj_list[i]) {
            head.push_back(to);
            weights.push_back(w);
        }
    }

    // Sentinel at the end
    first_out.push_back(head.size());
}

crp::Graph crp::Graph::reversed() const
{
    AdjacencyList adj(num_nodes());
    for (NodeId u = 0; u < num_nodes(); u++) {
        for (auto [v, weight] : (*this)[u]) {
            adj[v].push_back({u, weight});
        }
    }

    return Graph{adj};
}
