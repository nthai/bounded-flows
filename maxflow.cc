#include <boost/config.hpp>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <iomanip>

#include <boost/graph/edmonds_karp_max_flow.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/read_dimacs.hpp>
#include <boost/graph/graph_utility.hpp>

typedef boost::adjacency_list_traits<boost::vecS, boost::vecS, boost::directedS> Traits;
typedef boost::adjacency_list < boost::listS, boost::vecS, boost::directedS,
    boost::property<boost::vertex_name_t, std::string>,
    boost::property<boost::edge_capacity_t, double,
    boost::property<boost::edge_residual_capacity_t, double,
    boost::property<boost::edge_reverse_t, Traits::edge_descriptor>>>> Graph;

void print_graph(std::vector<std::vector<double>> const& g) {
    for (auto const& r : g) {
        for (auto const& x : r) {
            std::cout << std::setw(4) << x << " ";
        }
        std::cout << "\n";
    }
}

void print_graph(Graph& g) {
    boost::graph_traits<Graph>::vertex_iterator u_iter, u_end;
    boost::graph_traits<Graph>::out_edge_iterator ei, e_end;

    boost::property_map<Graph, boost::edge_capacity_t>::type capacity = boost::get(boost::edge_capacity, g);
    
    for (boost::tie(u_iter, u_end) = boost::vertices(g); u_iter != u_end; ++u_iter) {
        for (boost::tie(ei, e_end) = boost::out_edges(*u_iter, g); ei != e_end; ++ei) {
            std::cout << *u_iter << " " << boost::target(*ei, g) << " "
                << capacity[*ei] << std::endl;
        }
    }
}

std::vector<std::vector<double>> get_mod_capacities() {
    std::ifstream infile("graph.dat");
    int vcount, ecount;
    infile >> vcount >> ecount;

    std::vector<std::vector<double>> graph(vcount + 2, std::vector<double>(vcount + 2, 0));

    for (int i = 0; i < ecount; i++) {
        int from, to;
        double lower, upper;
        infile >> from >> to >> lower >> upper;

        graph[from + 1][to + 1] += upper - lower;
        graph[0][to + 1] += lower;
        graph[from + 1][vcount + 1] += lower;
    }

    // print_graph(graph);
    
    return graph;
}

int main() {
    std::vector<std::vector<double>> mod_capacities = get_mod_capacities();
    unsigned mod_vcount = mod_capacities.size();

    Graph g;

    std::vector<Traits::vertex_descriptor> vdesc(mod_vcount);
    for (int i = 0; i < mod_vcount; i++) {
        vdesc[i] = boost::add_vertex(g);
    }

    for (size_t from = 0; from < mod_vcount; from++) {
        for (size_t to = 0; to < mod_vcount; to++) {
            if (mod_capacities[from][to] != 0) {
                add_edge(from, to, { mod_capacities[from][to], 0.0 }, g);
            }
        }
    }
    // print_graph(g);

    Traits::vertex_descriptor s = vdesc[0], t = vdesc[mod_vcount - 1];
    boost::property_map<Graph, boost::edge_capacity_t>::type capacity = boost::get(boost::edge_capacity, g);
    boost::property_map<Graph, boost::edge_reverse_t>::type rev = boost::get(boost::edge_reverse, g);
    boost::property_map<Graph, boost::edge_residual_capacity_t>::type residual_capacity = boost::get(boost::edge_residual_capacity, g);
    std::vector<boost::default_color_type> color(boost::num_vertices(g));
    std::vector<Traits::edge_descriptor> pred(boost::num_vertices(g));
    double flow = boost::edmonds_karp_max_flow(g, s, t, capacity, residual_capacity, rev, &color[0], &pred[0]);
    std::cout << flow << std::endl;
}
