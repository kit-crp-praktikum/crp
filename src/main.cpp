#include "data-types.h"
#include "graph.h"
#include "lib/debug.h"
#include "lib/timer.h"

int main(int argc, char **argv)
{
    crp::Graph g(argv[1], argv[2], argv[3]);
    measure_time("Graph direct iteration", [=]() {
        int sumA = 0;
        int sumB = 0;

        for (int i = 0; i < g.num_nodes(); i++)
        {
            for (int j = g.first_out[i]; j < g.first_out[i + 1]; j++)
            {
                sumA += g.head[j] * i;
                sumB += g.weights[j] * i;
            }
        }

        std::cout << sumA _ sumB << std::endl;
    });

    measure_time("Iterator graph", [=]() {
        int sumA = 0;
        int sumB = 0;

        for (int i = 0; i < g.num_nodes(); i++)
        {
            for (auto [to, w] : g[i])
            { // edge i -> to with wieght w
                sumA += to * i;
                sumB += w * i;
            }
        }

        std::cout << sumA _ sumB << std::endl;
    });
}
