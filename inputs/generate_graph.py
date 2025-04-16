import random

def generate_clique(clique_size, min_capacity=1, max_capacity=10):
    vertices = list(range(clique_size))
    edges = []
    
    for i in vertices:
        for j in vertices:
            if (j == 1 and i != j) or (i != 1 and i < j):
                capacity = random.randint(min_capacity, max_capacity)
                edges.append((i, j, capacity))
    
    num_edges = len(edges)
    
    lines = [f"{clique_size} {num_edges}"]
    for edge in edges:
        lines.append(f"{edge[0]} {edge[1]} {edge[2]}")
    
    return "\n".join(lines)

def generate_sparse_graph(n, d, min_capacity=1, max_capacity=10):
    k = 2 # extra factor for source/sink

    if k * d > n - 1:
        raise ValueError("d must be at most n-1 (no self-loops allowed)")

    graph = {i: {} for i in range(n)}

    def has_reverse_edge(i, j):
        return i in graph[j]
    
    cap = random.randint(min_capacity, max_capacity)
    graph[0][1] = cap

    def add_edge(i, j):
        cap = random.randint(min_capacity, max_capacity)
        graph[i][j] = cap

    for i in range(n):
        if i == 1 :
            continue

        possible = set(range(n)) - {i} - {0}
        possible -= set(graph[i].keys())
        possible = {j for j in possible if not has_reverse_edge(i, j)}
        neighbors = random.sample(list(possible), (k * d if i == 0 else d))
        for v in neighbors:
            add_edge(i, v)

    possible = set(range(n)) - {0} - {1}
    possible = {j for j in possible if not has_reverse_edge(1, j)}
    extra_sample_num = k * d - (n - 2 - len(possible))
    if extra_sample_num > 0:
        neighbors = random.sample(list(possible), extra_sample_num)
        for v in neighbors:
            add_edge(v, 1)
    

    total_edges = sum(len(targets) for targets in graph.values())
    lines = [f"{n} {total_edges}"]
    for i in range(n):
        for j, cap in graph[i].items():
            lines.append(f"{i} {j} {cap}")
    return "\n".join(lines)

if __name__ == "__main__":
    num_vertices = 50000    
    degree = 500          
    min_cap = 1
    max_cap = 999

    try:
        graph_text = generate_sparse_graph(num_vertices, degree, min_cap, max_cap)
        with open(f"sparse/sparse_{num_vertices}_{degree}.txt", "w+") as f:
            f.write(graph_text)
    except ValueError as e:
        print("Error generating graph:", e)