#include <iostream>
#include <vector>
#include <queue>
#include <algorithm>
#include <climits>
#include <ctime>
#include <chrono>
#include <unistd.h>
#include <fstream>
#include <iomanip>
#include <string>
#include <atomic>

#include <omp.h>



using namespace std;

inline double to_sec(std::chrono::nanoseconds ns)
{
    return ns.count() * 1e-9;           // 1 ns  =  1e-9 s
}

#define SOURCE 0
#define SINK   1

static std::chrono::nanoseconds b_total{0};
static std::chrono::nanoseconds d_total{0};

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

inline int atomic_read(const int& x)
{
    int v;
    #pragma omp atomic read
    v = const_cast<int&>(x);
    return v;
}

inline int fwd_res(const Edge& e) { return e.capacity - atomic_read(e.flow); }
inline int rev_res(const Edge& e) { return atomic_read(e.flow); }


// Build level graph
bool construct_level_graph(Graph& graph, vector<int>& level) {
    int n = graph.size();
    level.assign(n, -1);
    level[SOURCE] = 0;
    queue<int> queue;
    queue.push(SOURCE);

    auto& adj = graph.getAdj();
    while (!queue.empty()) {
        int from = queue.front(); queue.pop();
        if (from == SINK) return true;
        for (const Edge& edge : adj[from])
            if (level[edge.to] < 0 && edge.flow < edge.capacity) {
                level[edge.to] = level[from] + 1;
                queue.push(edge.to);
            }
    }
    return level[SINK] >= 0;
}

bool construct_level_graph_parallel(Graph& g, std::vector<int>& level)
{
    const int n = g.size();
    level.assign(n, -1);
    level[SOURCE] = 0;

    std::vector<int> frontier{SOURCE}, next_frontier;
    auto& adj = g.getAdj();
    int d = 0;

    while (!frontier.empty())
    {
        ++d;
        next_frontier.clear();
        next_frontier.reserve(frontier.size() * 2);

        #pragma omp parallel
        {
            std::vector<int> local;
            local.reserve(256);

            #pragma omp for schedule(dynamic,1)
            for (size_t idx = 0; idx < frontier.size(); ++idx) {
                int u = frontier[idx];
                for (const Edge& e : adj[u]) {
                    if (fwd_res(e) == 0) continue;
                    int v = e.to;
                    if (level[v] != -1) continue;

                    //std::atomic_ref<int> lv(level[v]);
                    int expected = -1;
                    if (__sync_bool_compare_and_swap(&level[v], expected, d))
                        local.push_back(v);
                }
            }
            // merge thread-local buffer
            #pragma omp critical
            next_frontier.insert(next_frontier.end(),
                                 local.begin(), local.end());
        }

        if (level[SINK] != -1) return true;   
        frontier.swap(next_frontier);
    }
    return false;
}

// Send blocking flow (current‑arc optimisation)
int compute_block_flow(Graph& graph, int u, int pushed,
        vector<int>& level, vector<int>& next) {
    if (u == SINK || pushed == 0) return pushed;
    auto& adj = graph.getAdj();

    int through = 0;

    for (int& i = next[u]; i < (int)adj[u].size(); i++) {

        int j = i;
        if (j >= (int)adj[u].size()) break;

        Edge& edge = adj[u][j];
        if (edge.capacity > edge.flow && level[edge.to] == level[u] + 1) {
            int tr = compute_block_flow(graph, edge.to,
                         min(pushed - through, edge.capacity - edge.flow),
                         level, next);

            edge.flow += tr;
            adj[edge.to][edge.rev].flow -= tr;

            if (pushed == (through += tr)) {
                return through;
            }
        }
    }
    return through;
}

// main Dinic function
int maxFlow(Graph& graph) {
    if (SOURCE == SINK) return 0;
    int total = 0, n = graph.size();
    vector<int> level(n), next(n);

    double tb = 0.0;
    double td = 0.0;

    int loop_count = 0;
    while (true) {
        loop_count++;
        const auto tb0 = chrono::steady_clock::now();
        if (!(construct_level_graph(graph, level))) {
          break;
        }
        fill(next.begin(), next.end(), 0);
        const auto tb1 = chrono::steady_clock::now();
        tb += chrono::duration_cast<chrono::duration<double>>(tb1 - tb0).count();
        
        const auto td0 = chrono::steady_clock::now();
        while (int pushed = compute_block_flow(graph, SOURCE, INT_MAX, level, next))
            total += pushed;
        const auto td1 = chrono::steady_clock::now();
        td += chrono::duration_cast<chrono::duration<double>>(td1 - td0).count();
    }

    //cout << "Loops Count: " << loop_count << "\n";

    cout << "Phase 1 time: " << fixed << setprecision(3) << tb << '\n';
    cout << "Phase 2 time: " << fixed << setprecision(3) << td << '\n';

    return total;
}

int compute_block_flow_parallel_coarse(Graph& graph, int u, int pushed,
        vector<int>& level, vector<int>& next) {
    if (u == SINK || pushed == 0) return pushed;
    auto& adj = graph.getAdj();

    int end = (int)adj[u].size();

    int through = 0;

    int i;

    while(true) {
        #pragma omp atomic capture
        i = next[u]++;

        if (i >= end) break;

        Edge& edge = adj[u][i];
        if (edge.capacity > edge.flow && level[edge.to] == level[u] + 1) {
            int tr = compute_block_flow(graph, edge.to,
                         min(pushed - through, edge.capacity - edge.flow),
                         level, next);

            #pragma omp atomic
            edge.flow += tr;
            #pragma omp atomic
            adj[edge.to][edge.rev].flow -= tr;

            if (pushed == (through += tr)) {
                return through;
            }
        }
    }
    return through;
}

int compute_block_flow_parallel_start_coarse(Graph& graph, int u, int pushed,
        vector<int>& level, vector<int>& next, int num_threads) {
    if (u == SINK || pushed == 0) return pushed;
    auto& adj = graph.getAdj();

    const int start = next[u];
    const int end = (int)adj[u].size();

    int through = 0;

    #pragma omp parallel for schedule(dynamic, 4) num_threads(num_threads)
    for (int i = start; i < end; ++i) {
        Edge& edge = adj[u][i];
        if (edge.capacity > edge.flow && level[edge.to] == level[u] + 1) {
            int tr = compute_block_flow_parallel_coarse(graph, edge.to,
                         min(pushed - through, edge.capacity - edge.flow),
                         level, next);

            edge.flow += tr;
            adj[edge.to][edge.rev].flow -= tr;

            #pragma omp atomic
            through += tr;
        }
    }

    #pragma omp barrier
    
    return through;
}


int maxFlow_parallel_coarse(Graph& graph, int num_threads) {
    if (SOURCE == SINK) return 0;
    int total = 0, n = graph.size();
    vector<int> level(n), next(n);

    double tb = 0.0;
    double td = 0.0;

    int loop_count = 0;

    while (true) {
        loop_count++;

        const auto tb0 = chrono::steady_clock::now();
        if (!(construct_level_graph_parallel(graph, level))) {
          break;
        }
        fill(next.begin(), next.end(), 0);
        const auto tb1 = chrono::steady_clock::now();
        tb += chrono::duration_cast<chrono::duration<double>>(tb1 - tb0).count();
        
        const auto td0 = chrono::steady_clock::now();
        int break_counter = 0;
        while (true) {
            int pushed = compute_block_flow_parallel_start_coarse(graph, SOURCE, INT_MAX, level, next, num_threads);

            if (pushed == 0) {
                break_counter++;
                fill(next.begin(), next.end(), 0);
            } else {
                break_counter = 0;
            }

            total += pushed;

            if (break_counter >= 1) {
                break;
            } 
        }
            
        const auto td1 = chrono::steady_clock::now();
        td += chrono::duration_cast<chrono::duration<double>>(td1 - td0).count();
    }

    //cout << "Loops Count: " << loop_count << "\n";

    cout << "Phase 1 time: " << fixed << setprecision(3) << tb << '\n';
    cout << "Phase 2 time: " << fixed << setprecision(3) << td << '\n';

    return total;
}

inline bool reserve_push(Edge& fwd, Edge& rev, int delta)
{
    int old;
    #pragma omp atomic capture
    {
    old = fwd.flow;
    fwd.flow += delta;
    }

    if (old + delta > fwd.capacity) {
        #pragma omp atomic
        fwd.flow -= delta;
        return false;
    }
    #pragma omp atomic
    rev.flow -= delta;
    return true;
}

inline void refund(Edge& fwd, Edge& rev, int delta)
{
    #pragma omp atomic
    fwd.flow -= delta;
    #pragma omp atomic
    rev.flow += delta;
}


int flow_reserve(Graph& g,int u,int pushed,
                const std::vector<int>& level,std::vector<int>& next)
{
    if(u==SINK||pushed==0) return pushed;
    auto& adj = g.getAdj();
    for(int& i = next[u]; i<(int)adj[u].size(); ++i){
        Edge& fwd = adj[u][i];
        if(fwd_res(fwd)==0 || level[fwd.to]!=level[u]+1) continue; // saturated or not in level edge
        Edge& rev = adj[fwd.to][fwd.rev];

        int send = min(pushed, fwd_res(fwd));
        if(!reserve_push(fwd,rev,send)) continue;

        int got = flow_reserve(g,fwd.to,send,level,next);
        if(got==0){ refund(fwd,rev,send); continue; } // nothing got through
        if(got<send) refund(fwd,rev,send-got); // refund not needed flow
        return got;
    }
    return 0;
}

int maxFlow_parallel_fine(Graph& g)
{
    int n = g.size();
    vector<int> level(n);
    int total = 0;

    while (true) {
        auto t_b0 = chrono::steady_clock::now();
        bool reachable = construct_level_graph_parallel(g, level);
        b_total += chrono::steady_clock::now() - t_b0;
        if (! reachable) break;

        int phase_flow = 0;

        auto t_d0 = chrono::steady_clock::now();
        #pragma omp parallel reduction(+:phase_flow)
        {
            vector<int> next(n,0);

            int pushed;
            do {
                pushed = flow_reserve(g, SOURCE, INT_MAX, level, next);
                phase_flow += pushed;
            } while (pushed);
        }
        d_total += chrono::steady_clock::now() - t_d0;
        total += phase_flow;
    }
    return total;
}


int main(int argc, char* argv[]) {
    int num_threads = 1;
    std::string mode = "coarse"; 

    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];

        if (arg == "-n" && i + 1 < argc) {
            num_threads = std::stoi(argv[++i]);
        }
        else if (arg == "-m" && i + 1 < argc) {
            mode = argv[++i];
        }
    }

    omp_set_num_threads(num_threads);

    try {
        if (mode == "coarse") {
            const auto t0 = chrono::steady_clock::now();

            Graph g = loadGraph(cin);
            int result = 0;
            const auto t1 = chrono::steady_clock::now();
            const auto t2 = chrono::steady_clock::now();
            result = maxFlow_parallel_coarse(g, num_threads);
            const auto t3 = chrono::steady_clock::now();
            
            cout << "Max flow: " << result << '\n';
            const double load_time = chrono::duration_cast<chrono::duration<double>>(t1 - t0).count();
            const double compute_time = chrono::duration_cast<chrono::duration<double>>(t3 - t2).count();

            cout << "Loading time (sec): " << fixed << setprecision(3) << load_time << '\n';
            cout << "Computation time (sec): " << fixed << setprecision(3) << compute_time << '\n';
        }
        else {
            auto t0 = std::chrono::steady_clock::now();
            Graph g=loadGraph(cin);
            auto t1 = std::chrono::steady_clock::now();
            int flow=maxFlow_parallel_fine(g);
            auto t2 = std::chrono::steady_clock::now();

            cout << "Max flow: " << flow << '\n';
            const double load_time = chrono::duration_cast<chrono::duration<double>>(t1 - t0).count();
            const double compute_time = chrono::duration_cast<chrono::duration<double>>(t2 - t1).count();

            cout << "Phase 1 time: " << fixed << setprecision(3) << to_sec(b_total) << '\n';
            cout << "Phase 2 time: " << fixed << setprecision(3) << to_sec(d_total) << '\n';

            cout << "Loading time (sec): " << fixed << setprecision(3) << load_time << '\n';
            cout << "Computation time (sec): " << fixed << setprecision(3) << compute_time << '\n';

        }
    } catch (const exception& e) {
        cerr << "Error: " << e.what() << '\n';
        return 1;
    }
    return 0;
}
