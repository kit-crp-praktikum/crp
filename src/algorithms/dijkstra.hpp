#pragma once

#include <vector>
#include "src/datastructure/timestamped_vector.hpp"
#include "src/data-types.h"
#include "lib/id_queue.h"

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
class Dijkstra {
    public:
        Dijkstra(std::size_t size) : 
            visited(size, false),
            distance(size, INF),
            parent(size, INF),
            priority_queue(size),
            progress(0) {}

        //runs dijkstra until priority_queue is empty
        template<bool update_parents = false>
        void compute_distance(NodeId start, auto neighbors) {
            reset();
            init_start_node(start);
            while(!priority_queue.empty()) {
                step<update_parents>(neighbors);
            }
        }

        //runs dijkstra until priority_queue is empty or target is settled
        template<bool update_parents = false>
        void compute_distance_target(NodeId start, NodeId target, auto neighbors) {
            reset();
            init_start_node(start);
            NodeId settled = visited.size() + 2;
            while(!priority_queue.empty() && settled != target) {
                settled = step<update_parents>(neighbors);
            }
        }

        //caller must check that priority queue is not empty
        template<bool update_parents = false>
        NodeId step(auto neighbors) {
            assert(!priority_queue.empty());

            auto v = priority_queue.peek().id;
            auto dist = priority_queue.pop().key;

            visited[v] = true;
            distance[v] = dist;
            progress = dist;

            auto relax_operation = [&](NodeId u, Distance weight) {
                if(!visited[u]) {
                    Distance relaxed = distance[v] + weight;
                    distance[u] = std::min(distance[u], relaxed);
                    if(priority_queue.contains_id(u)) {
                        if constexpr(update_parents) {
                            if(relaxed < priority_queue.get_key(u)) {
                                parent[u] = v;
                            }
                        }
                        priority_queue.decrease_key({u, relaxed});  
                    } 
                    else {
                        if constexpr(update_parents) parent[u] = v;
                        priority_queue.push({u, relaxed});
                    }
                }
            };
            neighbors(v, relax_operation); //apply relax operation to all neighbors
            return v;
        }

        void reset() {
            visited.reset();
            distance.reset();
            parent.reset();
            priority_queue.clear();
            progress = 0;
        }

        void init_start_node(NodeId start) {
            priority_queue.push({start, 0});
        }

        bool priority_queue_is_empty() {
            return priority_queue.empty();
        }

        Distance get_progress() {
            return progress;
        }
        //Return distance to start node, if v was settled, otherwise INF.
        Distance tentative_distance(NodeId v) {
            return distance[v];
        }

        //Returns parent of node, if v was relaxed and computation was run with update_parents = true.
        //Otherwise contains INF. 
        NodeId get_parent(NodeId v) {
            return parent[v];
        }

    private:

        TimestampedVector<bool> visited;
        TimestampedVector<Distance> distance;
        TimestampedVector<NodeId> parent;
        MinIDQueue priority_queue;
        Distance progress; //last priority queue key settled
};

class BidirectionalDijstkra {
    public:
        BidirectionalDijstkra(std::size_t size) :
        number_of_nodes(size),
        fwd(size),
        bwd(size) {}

        //returns last settled node and distance between start and target
        template<bool update_parents = false>
        std::pair<NodeId, Distance> compute_distance_target(NodeId start, NodeId target, auto fwd_neighbors, auto bwd_neighbors) {
            fwd.reset(); bwd.reset();
            fwd.init_start_node(start); bwd.init_start_node(target);
            NodeId settled = number_of_nodes + 2;
            Distance tentative_distance = INF;
            std::size_t steps = 0;
            
            //alternating priority queue selection
            while(!(fwd.priority_queue_is_empty() && bwd.priority_queue_is_empty())
            && tentative_distance > fwd.get_progress() + bwd.get_progress()) {
                if((steps & 1 && !fwd.priority_queue_is_empty()) || bwd.priority_queue_is_empty()) {
                    settled = fwd.step<update_parents>(fwd_neighbors);
                }
                else {
                    settled = bwd.step<update_parents>(bwd_neighbors);
                }
                //no overlow, since INF = MAX_INT / 2
                tentative_distance = std::min(tentative_distance, fwd.tentative_distance(settled) + bwd.tentative_distance(settled));
                steps++;
            }
            return {settled, tentative_distance};
        }

        NodeId fwd_parent(NodeId v) {
            return fwd.get_parent(v);
        }

        NodeId bwd_parent(NodeId v) {
            return bwd.get_parent(v);
        }

    private:
        NodeId number_of_nodes;
        Dijkstra fwd;
        Dijkstra bwd;
};
