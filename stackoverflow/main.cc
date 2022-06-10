// based on: https://stackoverflow.com/a/67622091/7925855

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/boykov_kolmogorov_max_flow.hpp>
#include <boost/range/adaptors.hpp>
#include <fmt/ostream.h>
#include <fmt/ranges.h>
using boost::adaptors::filtered;

using Traits = boost::adjacency_list_traits<boost::vecS, boost::vecS, boost::directedS>;
using V      = Traits::vertex_descriptor;
using E      = Traits::edge_descriptor;

using Capacity = double;
using Color    = boost::default_color_type;

struct VertexProps {
    // std::string name;
    Color    color;
    Capacity distance;
    E        predecessor;
};

struct EdgeProps {
    int      id;
    Capacity weight, residual;
    E        reverse;
};

using Graph = boost::adjacency_list<
    boost::vecS, boost::vecS, boost::directedS,
    VertexProps,
    // see https://stackoverflow.com/a/64744086/85371 :(
    boost::property<boost::edge_capacity_t, Capacity, EdgeProps>>;

struct MyGraph {
    MyGraph(size_t nnodes) : _g(nnodes) {}

    void runSimulation(auto const& arcs, auto const& capacities)
    {
        reconfigure(arcs, capacities);

        Capacity maxflow = solve_max_flow(0, 5);

        auto cap       = get(boost::edge_capacity, _g);
        auto is_source = [this](V v) { return _g[v].color == Color::black_color; };

        fmt::print("Max flow {}\nNodes {} are in source subset\n", maxflow,
                vertices(_g) | filtered(is_source));

        for (E e : boost::make_iterator_range(edges(_g))) {
            bool st_cut =
                is_source(source(e, _g)) and
                not is_source(target(e, _g));

            fmt::print("Edge {} (id #{:2}), capacity {:3} {}\n", e, _g[e].id,
                    cap[e], st_cut ? "(ST Cut)" : "");
        }
    }

private:
    Graph _g;

    void reconfigure(auto const& arcs, auto const& capacities)
    {
        assert(arcs.size() == capacities.size());

        for (auto v : boost::make_iterator_range(vertices(_g))) {
            // boost::clear_out_edges(v, g);
            boost::clear_vertex(v, _g);
        }

        auto cap  = get(boost::edge_capacity, _g);
        auto eidx = get(&EdgeProps::id, _g);
        auto rev  = get(&EdgeProps::reverse, _g);

        auto  eindex = 0;

        for (auto [fr, to] : arcs) {
            auto edf  = add_edge(fr, to, _g).first;
            auto edr  = add_edge(to, fr, _g).first;
            eidx[edf] = 2 * eindex;
            eidx[edr] = eidx[edf] + 1;
            cap[edf]  = cap[edr] = capacities[eindex];

            rev[edf] = edr;
            rev[edr] = edf;

            ++eindex;
        }
    }

    Capacity solve_max_flow(V src, V sink)
    {
        return boykov_kolmogorov_max_flow(
            _g, src, sink,
            // named arguments
            boost::reverse_edge_map(get(&EdgeProps::reverse, _g))
                .residual_capacity_map(get(&EdgeProps::residual, _g))
                .vertex_color_map(get(&VertexProps::color,       _g))
                .predecessor_map(get(&VertexProps::predecessor,  _g))
                .distance_map(get(&VertexProps::distance,        _g))
            // end named arguments
        );
    }
};

int main() {
    MyGraph g{6};

    using namespace std;
    for (auto&& [arcs, capacities] : { tuple
            // 1
            {vector{pair{0, 1}, {0, 2}, {1, 2}, {1, 3}, {1, 4},
                                    {2, 4}, {3, 4}, {3, 5}, {4, 5}},
            vector{10, 10, 10, 10, 1, 4, 3, 2, 10}},
            // 2
            {vector{pair{0, 1}, {0, 2}, {1, 3}, {2, 4}, {3, 5}, {4, 5}},
            vector{10, 10, 10, 4, 2, 0}},
        })
    {
        g.runSimulation(arcs, capacities);
    }
}
