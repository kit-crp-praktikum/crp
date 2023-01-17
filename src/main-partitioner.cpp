#include <bitset>
#include <iostream>
#include <vector>
#include "data-types.h"
#include <iomanip>
#include <numeric>

// make private methods accessable for testing
#define private public
#include "crp/crp.h"
#include "tests/grid-graph.hpp"

int main()
{
    unsigned n = 32;
    std::vector<std::vector<std::pair<NodeId, Distance>>> gr(n);
    // path
    for (NodeId i = 0; i < n - 1; i++)
    {
        gr[i].push_back({i + 1, 1});
        gr[i + 1].push_back({i, 1});
    }
}
