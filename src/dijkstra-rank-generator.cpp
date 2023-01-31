#include "algorithms/dijkstra.hpp"
#include "cmdline-parsing.hpp"
#include "data-types.h"
#include "graph.h"
#include "lib/vector_io.h"
#include <cstdlib>
#include <filesystem>
#include <iostream>
#include <random>

#include <iostream>
#define _ << " " <<
#define debug(x) #x << " = " << x

int main(int argc, char **argv)
{
    int pos;

    pos = find_required_argument(argc, argv, 'i', "input", true);
    std::filesystem::path root_dir = argv[pos + 1];

    pos = find_required_argument(argc, argv, 'o', "output", true);
    std::filesystem::path output_dir = argv[pos + 1];

    pos = find_required_argument(argc, argv, 'w', "weight", true);
    std::string weight_type = argv[pos + 1];

    pos = find_required_argument(argc, argv, 'n', "number-starts", true);
    int number_of_starts = parse_integer_or_bail(argv[pos + 1]);

    check_input_directory_valid(root_dir.generic_string(), weight_type, false);

    crp::Graph g((root_dir / "first_out").generic_string(), (root_dir / "head").generic_string(), (root_dir / weight_type).generic_string());

    const int n = g.num_nodes();
    const int nr_ranks = (n <= 1 ? 0 : 31 - __builtin_clz(n - 1)) + 1;

    std::cout << debug(n) _ debug(nr_ranks) << std::endl;

    std::vector<std::vector<NodeId>> query_start(nr_ranks);
    std::vector<std::vector<NodeId>> query_end(nr_ranks);
    std::vector<std::vector<Distance>> answers(nr_ranks);

    std::random_device rd;
    std::mt19937 mt(rd());
    for (int i = 0; i < number_of_starts; i++)
    {
        NodeId start = std::uniform_int_distribution<NodeId>(0, n - 1)(mt);

        Dijkstra dijkstra(n);
        dijkstra.compute_distance(start, [&](NodeId u, auto relax_op) {
            for (auto [v, w] : g[u])
            {
                relax_op(v, w);
            }
        });

        std::vector<NodeId> nodes(n);
        std::iota(nodes.begin(), nodes.end(), 0u);
        std::sort(nodes.begin(), nodes.end(),
                  [&](NodeId a, NodeId b) { return dijkstra.tentative_distance(a) < dijkstra.tentative_distance(b); });

        for (int rank = 0; rank < nr_ranks; rank++)
        {
            NodeId end = nodes[(1 << rank)];

            std::cout << start _ end _ dijkstra.tentative_distance(end) << std::endl;

            query_start[rank].push_back(start);
            query_end[rank].push_back(end);
            answers[rank].push_back(dijkstra.tentative_distance(end));
        }
    }

    // Serialize in one array for storage
    std::vector<NodeId> start, end, ans;
    for (int i = 0; i < nr_ranks; i++)
    {
        for (auto &x : query_start[i])
            start.push_back(x);
        for (auto &x : query_end[i])
            end.push_back(x);
        for (auto &x : answers[i])
            ans.push_back(x);
    }

    save_vector((output_dir / "source").generic_string(), start);
    save_vector((output_dir / "target").generic_string(), end);
    save_vector((output_dir / (weight_type + "_length")).generic_string(), ans);
}
