#include <boost/config.hpp>
#include <fstream>

#include <boost/graph/edmonds_karp_max_flow.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graph_utility.hpp>

typedef boost::adjacency_list_traits<boost::vecS, boost::vecS, boost::directedS> Traits;
typedef boost::adjacency_list<boost::listS, boost::vecS, boost::directedS,
    boost::property<boost::vertex_name_t, std::string>,
    boost::property<boost::edge_capacity_t, long,
    boost::property<boost::edge_residual_capacity_t, long,
    boost::property<boost::edge_reverse_t, Traits::edge_descriptor>>>> Graph;

int main() {
    std::ifstream infile("network.dat");
    int vcount, ecount;
    infile >> vcount >> ecount;
    
    Graph g;
    std::vector<Traits::vertex_descriptor> vdesc(vcount);
    for (int i = 0; i < vcount; i++)
        vdesc[i] = boost::add_vertex(g);
    
    for (int i = 0; i < ecount; i++) {
        int from, to, capacity;
        infile >> from >> to >> capacity;
        boost::add_edge(vdesc[from], vdesc[to], { capacity, 0 }, g);
    }

    Traits::vertex_descriptor s = boost::vertex(0, g);
    Traits::vertex_descriptor t = boost::vertex(vcount - 1, g);

    boost::property_map<Graph, boost::edge_capacity_t>::type capacity = boost::get(boost::edge_capacity, g);
    boost::property_map<Graph, boost::edge_reverse_t>::type rev = boost::get (boost::edge_reverse, g);
    boost::property_map<Graph, boost::edge_residual_capacity_t>::type residual_capacity = get(boost::edge_residual_capacity, g);
    std::vector<boost::default_color_type> color(boost::num_vertices(g));
    std::vector<Traits::edge_descriptor> pred(boost::num_vertices(g));

    long flow = boost::edmonds_karp_max_flow(g, s, t, capacity, residual_capacity, rev, &color[0], &pred[0]);
    std::cout << flow << std::endl;
}