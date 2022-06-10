// based on: https://stackoverflow.com/a/67622091/7925855

#pragma once

#include <boost/graph/adjacency_list.hpp>
#include <boost/range/adaptors.hpp>

using Traits = boost::adjacency_list_traits<boost::vecS, boost::vecS, boost::directedS>;
using V = Traits::vertex_descriptor;
using E = Traits::edge_descriptor;

using Capacity = double;
using Color = boost::default_color_type;

struct VertexProps {
    Color color;
    Capacity distance;
    E predecessor;
};

struct EdgeProps {
    int id;
    Capacity weight, residual;
    E reverse;
};

using Graph = boost::adjacency_list<
    boost::vecS, boost::vecS, boost::directedS, VertexProps,
    boost::property<boost::edge_capacity_t, Capacity, EdgeProps>>;

class ModGraph {
    Graph _original, _modified;
    std::vector<std::pair<V, V>> arcs;
    std::vector<Capacity> capacities;
public:
    ModGraph(size_t nnodes, auto arcs, auto capacities) :
        _original(nnodes), _modified(nnodes + 2), arcs(arcs), capacities(capacities) {}
    
};
