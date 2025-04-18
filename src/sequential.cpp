#include <iostream>
#include <vector>
#include <queue>
#include <algorithm>
#include <climits>
#include <ctime>

using namespace std;

double elapsed_ms(clock_t start, clock_t end)
{
    return 1000.0 * (end - start) / CLOCKS_PER_SEC;
}

#define SOURCE 0
#define SINK   1

struct Edge {
    int to;       // endpoint
    int capacity; // capacity of the edge
    int flow;     // current flow
    int rev;      // index of reverse edge in adj[to]
};

class Graph {
    int N;                                  // #vertices
    int M;                                  // #edges
    vector<vector<Edge>> adj;               // adjacency list
public:
    explicit Graph(int n, int m) : N(n), M(m), adj(n) {}
    void addEdge(int u, int v, int cap) {
        Edge fwd{v, cap, 0, static_cast<int>(adj[v].size())};
        Edge rev{u, 0,   0, static_cast<int>(adj[u].size())};
        adj[u].push_back(fwd);
        adj[v].push_back(rev);
    }
    int size() const { return N; }
    vector<vector<Edge>>& getAdj() { return adj; }
};

Graph loadGraph(std::istream& in)
{
    int n, m;
    if (!(in >> n >> m)) throw std::runtime_error("incorrect input format");

    Graph graph(n, m);
    for (int i = 0; i < m; ++i) {
        int u, v, cap;
        in >> u >> v >> cap;
        graph.addEdge(u, v, cap);
    }
    return graph;
}


// Build level graph
bool constru_level_graph(Graph& graph, vector<int>& level) {
    int n = graph.size();
    level.assign(n, -1);
    level[SOURCE] = 0;
    queue<int> que;
    que.push(SOURCE);

    auto& adj = graph.getAdj();
    while (!que.empty()) {
        int from = que.front(); que.pop();
        for (const Edge& edge : adj[from])
            if (level[edge.to] < 0 && edge.flow < edge.capacity) {
                level[edge.to] = level[from] + 1;
                que.push(edge.to);
            }
    }
    return level[SINK] >= 0;
}

// Send blocking flow (current‑arc optimisation)
int compute_block_flow(Graph& graph, int u, int pushed,
        vector<int>& level, vector<int>& next) {
    if (u == SINK || pushed == 0) return pushed;
    auto& adj = graph.getAdj();

    for (int& i = next[u]; i < (int)adj[u].size(); ++i) {
        Edge& edge = adj[u][i];
        if (edge.capacity > edge.flow && level[edge.to] == level[u] + 1) {
            int tr = compute_block_flow(graph, edge.to,
                         min(pushed, edge.capacity - edge.flow),
                         level, next);
            if (tr) {
                edge.flow += tr;
                adj[edge.to][edge.rev].flow -= tr;
                return tr;
            }
        }
    }
    return 0;
}

// main Dinic function
int maxFlow(Graph& graph) {
    if (SOURCE == SINK) return 0;
    int total = 0, n = graph.size();
    vector<int> level(n), next(n);

    while (constru_level_graph(graph, level)) {
        fill(next.begin(), next.end(), 0);
        while (int pushed = compute_block_flow(graph, SOURCE, INT_MAX, level, next))
            total += pushed;
    }
    return total;
}


int main() {
    try {
        clock_t t0 = std::clock();

        Graph g = loadGraph(cin);
        clock_t t1 = std::clock();
        clock_t t2 = std::clock();
        cout << maxFlow(g) << '\n';
        clock_t t3 = std::clock();
        std::cout << "loading time = " << elapsed_ms(t0, t1) << " ms\n";
        std::cout << "compute time = " << elapsed_ms(t2, t3) << " ms\n";
    } catch (const exception& e) {
        cerr << "Error: " << e.what() << '\n';
        return 1;
    }
    return 0;
}
