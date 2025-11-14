import networkx as nx
import matplotlib.pyplot as plt
from itertools import combinations


class LayeredGraph:
    def __init__(self, vertices, edges):
        self.vertices = vertices
        self.edges = edges
        self.graph = self._create_graph()

    def _create_graph(self):
        G = nx.DiGraph()

        for layer_idx, layer_vertices in enumerate(self.vertices):
            for vertex_idx, vertex in enumerate(layer_vertices):
                node_id = (layer_idx, vertex_idx)
                G.add_node(node_id, layer=layer_idx, label=vertex)

        for layer_idx, layer_edges in enumerate(self.edges):
            if layer_idx + 1 < len(self.vertices):  # Ensure target layer exists
                for edge in layer_edges:
                    j, k = edge
                    if j < len(self.vertices[layer_idx]) and k < len(self.vertices[layer_idx + 1]):
                        source = (layer_idx, j)
                        target = (layer_idx + 1, k)
                        G.add_edge(source, target)

        return G

    def draw(self):
        pos = {}
        labels = {}
        for layer_idx, layer_vertices in enumerate(self.vertices):
            for vertex_idx, vertex in enumerate(layer_vertices):
                node_id = (layer_idx, vertex_idx)
                pos[node_id] = (vertex_idx, -layer_idx)
                labels[node_id] = vertex

        nx.draw(self.graph, pos, labels=labels, with_labels=True, node_size=2000, node_color='lightblue', font_size=10,
                font_weight='bold', arrows=True)
        plt.show()

    def is_isomorphic(self, other):
        self_labels = sorted([v for layer in self.vertices for v in layer])
        other_labels = sorted([v for layer in other.vertices for v in layer])

        if self_labels != other_labels:
            return False

        def extract_edges(graph):
            return sorted([(graph.vertices[u[0]][u[1]], graph.vertices[v[0]][v[1]]) for u, v in graph.graph.edges])

        return extract_edges(self) == extract_edges(other)

    def create_subgraph(self, starting_vertices):
        sub_nodes = set()
        for layer_idx, vertex_idx in starting_vertices:
            if layer_idx < len(self.vertices) and vertex_idx < len(self.vertices[layer_idx]):
                sub_nodes.add((layer_idx, vertex_idx))
                sub_nodes.update(nx.descendants(self.graph, (layer_idx, vertex_idx)))

        subgraph_vertices = [[] for _ in range(len(self.vertices))]
        for (layer_idx, vertex_idx) in sorted(sub_nodes):
            if vertex_idx < len(self.vertices[layer_idx]):
                subgraph_vertices[layer_idx].append(self.vertices[layer_idx][vertex_idx])

        subgraph_edges = [[] for _ in range(len(self.edges))]
        for u, v in self.graph.edges:
            if u in sub_nodes and v in sub_nodes:
                subgraph_edges[u[0]].append((subgraph_vertices[u[0]].index(self.vertices[u[0]][u[1]]),
                                             subgraph_vertices[v[0]].index(self.vertices[v[0]][v[1]])))

        return LayeredGraph(subgraph_vertices, subgraph_edges)

    def create_quotient_graph(self, target_vertices):
        quotient_nodes = set()
        for layer_idx, vertex_idx in target_vertices:
            if layer_idx < len(self.vertices) and vertex_idx < len(self.vertices[layer_idx]):
                quotient_nodes.add((layer_idx, vertex_idx))
                quotient_nodes.update(nx.ancestors(self.graph, (layer_idx, vertex_idx)))

        quotient_vertices = [[] for _ in range(len(self.vertices))]
        for (layer_idx, vertex_idx) in sorted(quotient_nodes):
            if vertex_idx < len(self.vertices[layer_idx]):
                quotient_vertices[layer_idx].append(self.vertices[layer_idx][vertex_idx])

        quotient_edges = [[] for _ in range(len(self.edges))]
        for u, v in self.graph.edges:
            if u in quotient_nodes and v in quotient_nodes:
                quotient_edges[u[0]].append((quotient_vertices[u[0]].index(self.vertices[u[0]][u[1]]),
                                             quotient_vertices[v[0]].index(self.vertices[v[0]][v[1]])))

        return LayeredGraph(quotient_vertices, quotient_edges)

    def generate_all_subgraphs(self):
        all_subgraphs = []
        all_vertices = [(layer_idx, v_idx) for layer_idx, layer in enumerate(self.vertices) for v_idx in
                        range(len(layer))]

        for subset in range(1, len(all_vertices) + 1):
            for sub_nodes in combinations(all_vertices, subset):
                subgraph = self.create_subgraph(sub_nodes)
                if not any(subgraph.is_isomorphic(existing) for existing in all_subgraphs):
                    all_subgraphs.append(subgraph)
        return all_subgraphs

    def generate_all_quotient_graphs(self):
        all_quotient_graphs = []
        all_vertices = [(layer_idx, v_idx) for layer_idx, layer in enumerate(self.vertices) for v_idx in
                        range(len(layer))]

        for subset in range(1, len(all_vertices) + 1):
            for target_nodes in combinations(all_vertices, subset):
                quotient_graph = self.create_quotient_graph(target_nodes)
                if not any(quotient_graph.is_isomorphic(existing) for existing in all_quotient_graphs):
                    all_quotient_graphs.append(quotient_graph)
        return all_quotient_graphs

    def generate_functions(self, other):
        functions = []
        all_quotients = self.generate_all_quotient_graphs()
        all_subgraphs = other.generate_all_subgraphs()

        for q in all_quotients:
            for s in all_subgraphs:
                if q.is_isomorphic(s):
                    functions.append((q, s))
        return functions

    def generate_meta_graph(self):
        subgraphs = self.generate_all_subgraphs()
        meta_vertices = subgraphs
        meta_edges = []

        for i, sub1 in enumerate(subgraphs):
            for j, sub2 in enumerate(subgraphs):
                if i != j and sub1.generate_functions(sub2):
                    meta_edges.append((i, j))

        return LayeredGraph(meta_vertices, meta_edges)

# Example usage
# Example usage
if __name__ == "__main__":
    graph3 = LayeredGraph([
        ['A', 'B', 'C', 'A'],
        ['B', 'C', 'A'],
        ['A', 'B'],
        ['C']
    ], [
        [(0, 0), (1, 0)],
        [(1, 0), (2, 0), (1, 1)],
        [(2, 0), (3, 0)]
    ])

    meta_graph = graph3.generate_meta_graph()
    meta_graph.draw()
