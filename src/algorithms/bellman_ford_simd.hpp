#pragma once

#include "data-types.h"
#include "datastructure/timestamped_vector.hpp"
#include "lib/id_queue.h"
#include <iostream>
#include <vector>
#include <array>
#include <emmintrin.h>
#include <immintrin.h>
#include <smmintrin.h>

#define SIMD_LEN 8

/**
 * Implementation of Bellman-Ford Algorithm using SIMD instructions.
 * Computes 8 one-to-all queries at once.
 *
 * Takes as lambda a function that accepts a void(NodeId v, auto f) function.
 *
 * Example:
 * auto neighbors = [&](NodeId v, auto f) {
 *      for(auto [u, weight] : gr[v]) {
 *          f(u, weight);
 *      }
 *  };
 *  auto myf = [&](NodeId u, Distance weight) {
 *      std::cout << u << " " << weight << "\n";
 *  };
 *  neighbors(0, myf);
 *
 *
 * auto for_all_nodes = [&](auto f) {
 *     for(NodeId v = 0; v < n; v++) {
 *          f(v);
 *     }
 * };
 */

class BellmanFordSIMD
{
  public:
    BellmanFordSIMD(std::size_t size) : distance(size), number_of_nodes(size)
    {
    }

    void generic_compute_distance(std::array<NodeId, SIMD_LEN> &start_nodes, auto neighbors, auto for_all_nodes)
    {
        reset();
        init_start_nodes(start_nodes);
        bool changed = true;
        for (uint32_t i = 0; i < number_of_nodes - 1 && changed; i++)
        {
            changed = false;
            auto relax_neighbors = [&](NodeId v) {
                // TODO test aligned loading, different loading
                __m256i vec_v = _mm256_loadu_si256((__m256i_u *)(distance[v].data()));
                auto relax_operation = [&](NodeId u, Distance w) {
                    std::array<Distance, SIMD_LEN> e_w = {w, w, w, w, w, w, w, w};
                    // 8x 32-bit integer into simd register
                    __m256i vec_u = _mm256_loadu_si256((__m256i_u *)(distance[u].data()));
                    __m256i vec_weight = _mm256_loadu_si256((__m256i_u *)(e_w.data()));

                    // relax edge
                    __m256i vec_relaxed = _mm256_add_epi32(vec_v, vec_weight);
                    
                    // perform componentwise minimum for unsigned integers
                    __m256i vec_min = _mm256_min_epu32(vec_u, vec_relaxed);
                    _mm256_storeu_si256((__m256i_u *)distance[u].data(), vec_min);
                    
                    // check if atleast one entry changed
                    __m256i c = _mm256_cmpeq_epi32(vec_u, vec_min);
                    uint32_t mask = _mm256_movemask_epi8(c);
                    if (mask != 0xffffffff)
                        changed = true;
                };
                neighbors(v, relax_operation);
            };
            for_all_nodes(relax_neighbors);
        }
    }

    void compute_distance(std::array<NodeId, SIMD_LEN> &start_nodes, auto neighbors)
    {
        auto for_all_nodes = [&](auto f) {
            for (NodeId v = 0; v < number_of_nodes; v++)
            {
                f(v);
            }
        };
        generic_compute_distance(start_nodes, neighbors, for_all_nodes);
    }


    std::array<Distance, SIMD_LEN> tentative_distance(NodeId to)
    {
        return distance[to];
    }

    // set num nodes if you want to operate on a subgraph of smaller size. NodeIDs have to be mapped to 0, ..., n - 1
    void set_number_of_nodes(std::size_t n)
    {
        number_of_nodes = n;
    }

  private:
    void init_start_nodes(std::array<NodeId, SIMD_LEN> &start_nodes)
    {
        for (int i = 0; i < SIMD_LEN; i++)
        {
            distance[start_nodes[i]][i] = 0;
        }
    }

    void reset()
    {
        for (uint32_t i = 0; i < number_of_nodes; i++)
        {   
            std::fill(distance[i].data(), distance[i].data() + SIMD_LEN, INF);
        }
    }

    std::vector<std::array<Distance, SIMD_LEN>> distance;
    std::size_t number_of_nodes;
};
