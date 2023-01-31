#pragma once

#include "lib/id_queue.h"
#include "src/data-types.h"
#include "src/datastructure/timestamped_vector.hpp"
#include <vector>
#include <iostream>


/**
 * Implementation of Dijkstras Algorithm.
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
 */
class Dijkstra
{
  public:
    Dijkstra(std::size_t size) : distance(size, INF), parent(size, INF), priority_queue(size), progress(0)
    {
    }

    // runs dijkstra until priority_queue is empty
    template <bool update_parents = false> void compute_distance(NodeId start, auto neighbors)
    {
        reset();
        init_start_node(start);
        while (!priority_queue.empty())
        {
            step<update_parents>(neighbors);
        }
    }

    // runs dijkstra until priority_queue is empty or target is settled
    template <bool update_parents = false> void compute_distance_target(NodeId start, NodeId target, auto neighbors)
    {
        reset();
        init_start_node(start);
        NodeId settled = distance.size() + 2;
        while (!priority_queue.empty() && settled != target)
        {
            settled = step<update_parents>(neighbors);
        }
    }

    // caller must check that priority queue is not empty
    template <bool update_parents = false> NodeId step(auto neighbors)
    {
        assert(!priority_queue.empty());

        auto v = priority_queue.peek().id;
        auto dist = priority_queue.pop().key;

        distance[v] = dist;
        progress = dist;

        auto relax_operation = [&](NodeId u, Distance weight) {
            Distance relaxed = distance[v] + weight;
            if (relaxed < distance[u])
            {
                distance[u] = relaxed;
                if (priority_queue.contains_id(u))
                {
                    if constexpr (update_parents)
                    {
                        if (relaxed < priority_queue.get_key(u))
                        {
                            parent[u] = v;
                        }
                    }
                    #ifdef PRINT_EDGES
                        if (relaxed < priority_queue.get_key(u))
                        {
                            logging.push_back(u);
                            logging.push_back(v);
                        }
                    #endif
                    priority_queue.decrease_key({u, relaxed});
                }
                else
                {
                    if constexpr (update_parents)
                        parent[u] = v;
                    
                    #ifdef PRINT_EDGES
                        logging.push_back(u);
                        logging.push_back(v);
                    #endif
                    priority_queue.push({u, relaxed});
                }
            }
        };
        neighbors(v, relax_operation); // apply relax operation to all neighbors
        return v;
    }

    void reset()
    {
        distance.reset();
        parent.reset();
        priority_queue.clear();
        progress = 0;
    }

    void init_start_node(NodeId start)
    {
        priority_queue.push({start, 0});
    }

    bool priority_queue_is_empty()
    {
        return priority_queue.empty();
    }

    Distance get_progress()
    {
        return progress;
    }
    // Return distance to start node, if v was settled, otherwise INF.
    Distance tentative_distance(NodeId v)
    {
        return distance[v];
    }

    // Returns parent of node, if v was relaxed and computation was run with update_parents = true.
    // Otherwise contains INF.
    NodeId get_parent(NodeId v)
    {
        return parent[v];
    }

    std::vector<NodeId> unpack(NodeId s, NodeId t)
    {
        NodeId node = t;
        std::vector<NodeId> path;
        while (node != s)
        {
            path.push_back(node);
            node = get_parent(node);
        }
        path.push_back(node);
        std::reverse(path.begin(), path.end());
        return path;
    }

    // to print edges of dijkstra run
    #ifdef PRINT_EDGES
    std::vector<NodeId> logging;
    #endif
  private:
    TimestampedVector<Distance> distance;
    TimestampedVector<NodeId> parent;
    MinIDQueue priority_queue;
    Distance progress; // last priority queue key settled
};

class BidirectionalDijstkra
{
  public:
    BidirectionalDijstkra(std::size_t size) : number_of_nodes(size), fwd(size), bwd(size)
    {
    }

    // returns meeting point and distance between start and target
    template <bool update_parents = false>
    std::pair<NodeId, Distance> compute_distance_target(NodeId start, NodeId target, auto fwd_neighbors,
                                                        auto bwd_neighbors)
    {
        fwd.reset();
        bwd.reset();
        fwd.init_start_node(start);
        bwd.init_start_node(target);
        NodeId settled = number_of_nodes + 2;
        NodeId meeting_point = target;
        Distance tentative_distance = INF;
        std::size_t steps = 0;

        // alternating priority queue selection
        while (!(fwd.priority_queue_is_empty() && bwd.priority_queue_is_empty()) &&
               tentative_distance > fwd.get_progress() + bwd.get_progress())
        {
            if ((steps & 1 && !fwd.priority_queue_is_empty()) || bwd.priority_queue_is_empty())
            {
                settled = fwd.step<update_parents>(fwd_neighbors);
            }
            else
            {
                settled = bwd.step<update_parents>(bwd_neighbors);
            }
            // no overlow, since INF = MAX_INT / 2
            if (tentative_distance > fwd.tentative_distance(settled) + bwd.tentative_distance(settled))
            {
                tentative_distance = fwd.tentative_distance(settled) + bwd.tentative_distance(settled);
                meeting_point = settled;
            }
            steps++;
        }
        #ifdef PRINT_EDGES
            std::cerr << "fwd search: " << fwd.logging.size() << ", bwd search:" << bwd.logging.size() << "\n";
            std::vector<NodeId> log;
            log.push_back(fwd.logging.size());
            log.push_back(bwd.logging.size());
            log.insert(log.end(), fwd.logging.begin(), fwd.logging.end());
            log.insert(log.end(), bwd.logging.begin(), bwd.logging.end());
            fwd.logging.clear();
            bwd.logging.clear();
            std::cout.write((char *)log.data(), sizeof(uint32_t) * log.size());
        #endif
        
        return {meeting_point, tentative_distance};
    }

    NodeId fwd_parent(NodeId v)
    {
        return fwd.get_parent(v);
    }

    NodeId bwd_parent(NodeId v)
    {
        return bwd.get_parent(v);
    }

    std::vector<NodeId> unpack(NodeId s, NodeId t, NodeId meeting_point)
    {
        std::vector<NodeId> path;
        NodeId node = meeting_point;
        while (node != s)
        {
            path.push_back(node);
            node = fwd_parent(node);
        }
        path.push_back(node);
        std::reverse(path.begin(), path.end());
        // now fwd path added: s..meeting_point
        // meeting_point..t is still to be added
        node = meeting_point;
        while (node != t)
        {
            node = bwd_parent(node);
            path.push_back(node);
        }
        return path;
    }

  private:
    NodeId number_of_nodes;
    Dijkstra fwd;
    Dijkstra bwd;
};
