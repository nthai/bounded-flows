#include <iostream>
#include <vector>
#include <iomanip>
#include <fstream>
#include <queue>

using OriginalAdjMatrix = std::vector<std::vector<std::pair<double, double>>>;
using AdjMatrix = std::vector<std::vector<double>>;

OriginalAdjMatrix get_original_adj_matrix() {
    using namespace std;
    ifstream infile("../graph.dat");
    int vcount, ecount;
    infile >> vcount >> ecount;

    OriginalAdjMatrix graph(vcount, vector<pair<double, double>>(vcount, make_pair(0, 0)));

    for (int i = 0; i < ecount; i++) {
        int from, to;
        double lower, upper;
        infile >> from >> to >> lower >> upper;

        graph[from][to].first = lower;
        graph[from][to].second = upper;
    }

    return graph;
}

AdjMatrix get_modified_adj_matrix(OriginalAdjMatrix const& original) {
    using namespace std;
    
    size_t newsize = original.size() + 2;
    AdjMatrix newgraph = vector<vector<double>>(newsize, vector<double>(newsize, 0));

    for (size_t i = 0; i < original.size(); ++i) {
        for (size_t j = 0; j < original[i].size(); j++) {
            newgraph[0][j + 1] += original[i][j].first;
            newgraph[i + 1][newsize - 1] += original[i][j].first;
            newgraph[i + 1][j + 1] = original[i][j].second - original[i][j].first;
        }
    }

    newgraph[newsize - 2][1] = 9999;
    
    return newgraph;
}

void print_graph(AdjMatrix const& g) {
    for (auto& row : g) {
        for (auto& x : row) {
            std::cout << std::setw(4) << x << " ";
        }
        std::cout << "\n";
    }
}

bool bfs(AdjMatrix g, int s, int t, std::vector<int>& parent) {
    // std::cout << "Looking for path between " << s << " and " << t;
    
    bool visited[g.size()] = { 0 };
    std::queue<int> q;
    q.push(s);
    visited[s] = true;
    parent[s] = -1;
 
    while (!q.empty()) {
        int u = q.front();
        q.pop();
        for (int v = 0; v < g.size(); v++) {
            if (visited[v] == false && g[u][v] > 0) {
                q.push(v);
                parent[v] = u;
                visited[v] = true;
            }
        }
    }
    
    // std::cout << (visited[t] == true ? ": found\n" : ": not found\n");

    return (visited[t] == true);
}

std::tuple<double, AdjMatrix, AdjMatrix> maxflow(AdjMatrix g, int s, int t) {
    int u, v;
    AdjMatrix rGraph(g.size(), std::vector<double>(g.size(), 0));  
    AdjMatrix fGraph(g.size(), std::vector<double>(g.size(), 0));  
    for (u = 0; u < g.size(); u++) {
        for (v = 0; v < g.size(); v++) {
            rGraph[u][v] = g[u][v];
        }
    }
    std::vector<int> parent(g.size(), -1);
    double max_flow = 0;
    while (bfs(rGraph, s, t, parent)) {
        double path_flow = 99999;
        for (v = t; v != s; v = parent[v])
        {
            u = parent[v];
            path_flow = std::min(path_flow, rGraph[u][v]);
            // std::cout << "path flow: " << path_flow << std::endl;
        }
        for (v = t; v != s; v = parent[v])
        {
            u = parent[v];
            rGraph[u][v] -= path_flow;
            rGraph[v][u] += path_flow;

            fGraph[u][v] += path_flow;
            fGraph[v][u] -= path_flow;
        }
        max_flow += path_flow;
    }
    return { max_flow, rGraph, fGraph };
}

std::tuple<double, AdjMatrix, AdjMatrix> improve_flow(AdjMatrix residual, AdjMatrix fGraph, int s, int t) {
    double plus_flow = 0;
    std::vector<int> parent(residual.size(), -1);
    int u, v;
    while (bfs(residual, s, t, parent)) {
        std::cout << "FOUND!\n";
        double path_flow = 99999;
        for (v = t; v != s; v = parent[v])
        {
            u = parent[v];
            path_flow = std::min(path_flow, residual[u][v]);
            // std::cout << "path flow: " << path_flow << std::endl;
        }
        for (v = t; v != s; v = parent[v])
        {
            u = parent[v];
            residual[u][v] -= path_flow;
            residual[v][u] += path_flow;

            fGraph[u][v] += path_flow;
            fGraph[v][u] -= path_flow;
        }
        plus_flow += path_flow;
    }
    return { plus_flow, residual, fGraph };
}

void print_graph(OriginalAdjMatrix const& g) {
    for (auto& row : g) {
        for (auto& x : row) {
            std::cout << "<"
                << std::setw(4) << x.first << ", "
                << std::setw(4) << x.second << "> ";
        }
        std::cout << "\n";
    }
}

OriginalAdjMatrix set_flows(OriginalAdjMatrix const& orig, AdjMatrix const& flow) {
    OriginalAdjMatrix flows(orig);
    for (int i = 0; i < orig.size(); ++i) {
        for (int j = 0; j < orig.size(); ++j) {
            if (flow[i + 1][j + 1] > 0)
                flows[i][j].first += flow[i + 1][j + 1];
        }
    }
    return flows;
}

int main() {
    OriginalAdjMatrix orig = get_original_adj_matrix();
    AdjMatrix mat = get_modified_adj_matrix(orig);

    print_graph(mat);
    auto [maxf, rGraph, fGraph] = maxflow(mat, 0, 7);

    std::cout << maxf << std::endl;
    std::cout << "Residual graph: \n";
    print_graph(rGraph);
    print_graph(orig);
    
    std::cout << "New graph\n";
    auto flows = set_flows(orig, fGraph);
    print_graph(flows);

    double improved_flow;
    std::tie(improved_flow, rGraph, fGraph) = improve_flow(rGraph, fGraph, 0, 5);
    std::cout << "Residual\n";
    print_graph(rGraph);
}