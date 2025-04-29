import random

def generate_nn_graph(num_layers, nodes_per_layer, min_capacity=1, max_capacity=10):
    total_vertices = 2 + num_layers * nodes_per_layer

    def node_id(layer, pos):
        return 2 + (layer - 1) * nodes_per_layer + pos

    edges = []
    for j in range(nodes_per_layer):
        cap = random.randint(min_capacity, max_capacity)
        edges.append((0, node_id(1, j), cap))

    for layer in range(1, num_layers):
        for u in range(nodes_per_layer):
            for v in range(nodes_per_layer):
                cap = random.randint(min_capacity, max_capacity)
                edges.append((node_id(layer, u),
                              node_id(layer + 1, v),
                              cap))

    for j in range(nodes_per_layer):
        cap = random.randint(min_capacity, max_capacity)
        edges.append((node_id(num_layers, j), 1, cap))

    total_edges = len(edges)
    lines = [f"{total_vertices} {total_edges}"]
    for u, v, c in edges:
        lines.append(f"{u} {v} {c}")
    return "\n".join(lines)

if __name__ == "__main__":        
    min_cap = 1
    max_cap = 999

    num_layers = 20
    nodes_per_layer = 512

    try:
        graph_text = generate_nn_graph(num_layers, nodes_per_layer, min_cap, max_cap)
        with open(f"dense/nn_{nodes_per_layer}_{num_layers}.txt", "w+") as f:
            f.write(graph_text)

    except ValueError as e:
        print("Error generating graph:", e)