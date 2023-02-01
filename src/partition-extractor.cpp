#include "cmdline-parsing.hpp"
#include <lib/vector_io.h>

int main(int argc, char **argv)
{
    int pos;
    pos = find_required_argument(argc, argv, 'c', "cells-per-level", true);
    const int cells_per_level = parse_integer_or_bail(argv[pos + 1]);

    pos = find_required_argument(argc, argv, 'o', "output-levels", true);
    const int output_levels = parse_integer_or_bail(argv[pos + 1]);

    const int bits_per_level = 32 - __builtin_clz(cells_per_level - 1);
    const int total_bits = output_levels * bits_per_level;
    const uint32_t mask = (1u << total_bits) - 1;

    pos = find_required_argument(argc, argv, 'i', "input", true);

    auto partition = load_vector<uint32_t>(argv[pos + 1]);
    for (auto &x : partition)
        x &= mask;
    std::cout.write((char *)partition.data(), sizeof(uint32_t) * partition.size());
}
