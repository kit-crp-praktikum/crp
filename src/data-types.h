#pragma once

#include <cstdint>
#include <limits>
#include <vector>

using NodeId = uint32_t;
using Distance = uint32_t;
using Path = std::vector<NodeId>;

const unsigned invalid_id = 4294967295u;                       // for priority queue
const Distance INF = std::numeric_limits<Distance>::max() / 2; // to allow a + b without overflow
