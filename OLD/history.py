from itertools import combinations, chain
from venv import create

import networkx as nx
import matplotlib.pyplot as plt
from collections import Counter

from networkx.algorithms.isomorphism.isomorph import is_isomorphic

from graph_operations import create_radical_layers


#from itertools import combinations

#from networkx.algorithms.shortest_paths.unweighted import predecessor
#from networkx.classes import all_neighbors


def all_sublists(n: int):
    L = list(range(n))

    return [list(subset) for subset in chain.from_iterable(combinations(L, r) for r in range(1,len(L) + 1))]


class ModuleDiagram:
    def __init__(self, vertex_labels, arrows):
        self.vertex_labels = vertex_labels # list of strings; one for each vertex
        self.num_vertices = len(vertex_labels)
        self.arrows = arrows # list of tuples of vertex indices; (source,sink)
        self.num_arrows = len(arrows)
        self.basic_graph = self._create_basic_graph()
        self.radical_layers = self.create_radical_layers()
        self.radical_labels = self.create_radical_labels()
        self.graph = self._create_layered_graph()
        self.top = self.find_top()
        self.socle = self.find_socle()
        self.all_sub = self.generate_all_sub_modules()
        self.all_quot = self.generate_all_quotient_modules()

    def _create_basic_graph(self):
        G = nx.DiGraph()
        for vertex_id, vertex_label in enumerate(self.vertex_labels):
            G.add_node(vertex_id, label = vertex_label)
        G.add_edges_from(self.arrows)
        return G

    def create_radical_layers(self):  # G is a graph in networkx
        # A vertex is in radical layer l if it is l steps away from a vertex that is a source
        G = self.basic_graph
        sources = [x for x in G.nodes() if G.in_degree(x) == 0]
        if sources == []:
            print('This nx.DiGraph has no vertices that are sources.')
            return False
        radical_layers = [sources]
        test_vertices = sources
        vertices_left = list(G.nodes())
        while vertices_left != []:  # This is to make sure each vertex is assigned a unique radical level
            next_neighbors = []
            for vertex in test_vertices:
                next_neighbors += [x for x in G.neighbors(vertex) if x in vertices_left]
                vertices_left.remove(vertex)
            if next_neighbors != []:
                radical_layers += [next_neighbors]
                test_vertices = next_neighbors
            else:
                return radical_layers

    def create_radical_labels(self):
        radical_labels = []
        for layer in self.radical_layers:
            radical_labels += [[self.vertex_labels[i] for i in layer]]
        return radical_labels


    def find_top(self):
        # The top of a DiGraph is the set of vertices that are sources
        return self.create_radical_layers()[0]

    def find_socle(self):
        # The socle of a DiGraph is the set of vertices that are sinks
        G = self.basic_graph
        socle = [x for x in G.nodes() if G.out_degree(x) == 0]
        if socle == []:
            print('This nx.DiGraph has no vertices that are sinks.')
            return False
        else:
            return socle

    def relabel_mapping_wrt_radical_layers(self):
        # This changes the labelling of the vertices in the DiGraph so that the vertices in radical layer l are strictly smaller than those in radical layer l+1
        radical_layers = self.create_radical_layers()
        mapping_canonical_to_example = {}
        mapping_example_to_canonical = {}
        counter = 0
        for layer in radical_layers:
            for vertex in layer:
                mapping_canonical_to_example[counter] = vertex
                mapping_example_to_canonical[vertex] = counter
                counter += 1
        new_vertex_labels = [self.vertex_labels[mapping_canonical_to_example[i]] for i, j in
                             enumerate(self.vertex_labels)]
        new_arrows = [(mapping_example_to_canonical[i], mapping_example_to_canonical[j]) for (i, j) in self.arrows]
        return new_vertex_labels, new_arrows


    def _create_layered_graph(self): # takes the ModuleDiagram and interprets it as a networkx graph
        G = nx.DiGraph()
        for layer_idx, layer_vertices in enumerate(self.radical_layers):
            for vertex in layer_vertices:
                G.add_node(vertex, layer=layer_idx, label=(vertex,self.vertex_labels[vertex]))
        for arrow in self.arrows:
            G.add_edge(*arrow)
        return G

    def draw(self): # draws the networkx graph nicely
        pos = {}
        labels = {}
        for layer_idx, layer_vertices in enumerate(self.radical_layers):
            for vertex_idx, vertex in enumerate(layer_vertices):
                pos[vertex] = (vertex_idx, -layer_idx)
                labels[vertex] = vertex,self.vertex_labels[vertex]

        nx.draw(self.graph, pos, labels=labels, with_labels=True, node_size=1000, node_color='lightblue', font_size=10,
                font_weight='bold', arrows=True)
        plt.show()


    def create_sub_module(self, generating_vertices):
        if generating_vertices == []:
            return
        generating_vertices = list(set(generating_vertices))
        test_vertices = generating_vertices
        all_vertices = set(generating_vertices)
        while test_vertices != []:
            new_vertices = []
            for vertex in test_vertices:
                new_vertices += list(self.graph.successors(vertex))
                all_vertices.add(vertex)
            test_vertices = new_vertices
        all_vertices = list(all_vertices)
        mapping = {}
        for i, j in enumerate(all_vertices):
            mapping[j] = i
        inclusion_map = [(mapping[j], j) for j in all_vertices]
        H = nx.induced_subgraph(self.graph, all_vertices)
        H = nx.relabel_nodes(H, mapping)
        vertex_labels = [self.vertex_labels[i] for i in all_vertices]
        arrows = [a for a in H.edges]
        M = ModuleDiagram(vertex_labels, arrows)
        return [self, M, inclusion_map]

    def create_quotient_module(self, generating_vertices):
        if generating_vertices == []:
            return
        generating_vertices = list(set(generating_vertices))
        test_vertices = generating_vertices
        all_vertices = set(generating_vertices)
        while test_vertices != []:
            new_vertices = []
            for vertex in test_vertices:
                new_vertices += list(self.graph.predecessors(vertex))
                all_vertices.add(vertex)
            test_vertices = new_vertices
        all_vertices = list(all_vertices)
        mapping = {}
        for i, j in enumerate(all_vertices):
            mapping[j] = i
        quotient_map = [(mapping[j], j) for j in all_vertices]
        H = nx.induced_subgraph(self.graph, all_vertices)
        H = nx.relabel_nodes(H, mapping)
        vertex_labels = [self.vertex_labels[i] for i in all_vertices]
        arrows = [a for a in H.edges]
        M = ModuleDiagram(vertex_labels, arrows)
        return [self, M, quotient_map]

    def is_isomorphic_to(self,other):
        def node_match(n1, n2):
            return n1['layer'] == n2['layer'] and n1['label'][1] == n2['label'][1]
        return nx.is_isomorphic(self.graph, other.graph, node_match = node_match)

    def generate_all_sub_modules(self):
        all_submodules = [self.create_sub_module([0])]
        all_generating_vertices = all_sublists(self.num_vertices)
        for generating_vertices in all_generating_vertices:
            new_submodule = self.create_sub_module(generating_vertices)
            H = new_submodule[1].graph
            G = H.to_undirected()
            if nx.is_connected(G):
                add = True
                while add:
                    for submodule in all_submodules:
                        if new_submodule[1].is_isomorphic_to(submodule[1]):
                            add = False
                    if add:
                        all_submodules += [new_submodule]
        return all_submodules

    def generate_all_quotient_modules(self):
        all_quotient_modules = [self.create_quotient_module([0])]
        all_generating_vertices = all_sublists(self.num_vertices)
        for generating_vertices in all_generating_vertices:
            new_quotient_module = self.create_quotient_module(generating_vertices)
            H = new_quotient_module[1].graph
            G = H.to_undirected()
            if nx.is_connected(G):
                add = True
                while add:
                    for submodule in all_quotient_modules:
                        if new_quotient_module[1].is_isomorphic_to(submodule[1]):
                            add = False
                    if add:
                        all_quotient_modules += [new_quotient_module]
        return all_quotient_modules

    def homomorphism_group(self,other):
        hom = []
        for quotient in self.all_quot:
            print(quotient)
            quotient = quotient[1]
            for sub in other.all_sub:
                sub = sub[1]
                if quotient.is_isomorphic_to(sub):
                    hom += [self,quotient,sub,other]
        return hom


class SubModuleDiagram:

    def __init__(self, super_module, sub_module, inclusion_map):
        self.super_module = super_module # a ModuleDiagram object
        self.sub_module = sub_module # a ModuleDiagram object
        self.inclusion_map = inclusion_map # a list of tuples




# Example usage
if __name__ == "__main__": # can import it to another file and nothing below this line will run
    # this is where you put test cases so you can see what the code does

    M = ModuleDiagram(['1','1','3','2'], [(0,1),(0,2),(2,3)])
    L = ModuleDiagram(['2','1','1','3'], [(2,1),(2,3),(3,0)])
    N = ModuleDiagram(['2', 'A', '1', '3'], [(2, 1), (2, 3), (3, 0)])

    #M.draw()
    #print(M.create_sub_module([3]))
    #print('hi',M.create_sub_module([3])[1].radical_labels)
    #print(M.create_quotient_module([2,2]))
    #print(M.top)
    #print(M.socle)
    #print('nodes ', M.graph.nodes[2]['layer'], M.graph.nodes[2]['label'][1])
    #print('isomorphism ', M.is_isomorphic_to(N))
    #print('isomorphism ', M.is_isomorphic_to(L))
    #print('submodules', len(M.generate_all_sub_modules()))
    print('submodules', [M.generate_all_sub_modules()[i][2] for i in range(len(M.generate_all_sub_modules()))])
    print('quotient modules', [M.generate_all_quotient_modules()[i][2] for i in range(len(M.generate_all_quotient_modules()))])
    print('original vertex labels ', M.vertex_labels)
    print('radical layers ', M.radical_layers)
    print('radical labels', M.radical_labels)
    #print('hom', M.homomorphism_group(M))






#
#
#
#     def is_isomorphic(self, other):
#         self_labels = sorted([v for layer in self.vertices for v in layer])
#         other_labels = sorted([v for layer in other.vertices for v in layer])
#
#         if self_labels != other_labels:
#             return False
#
#         def extract_edges(graph):
#             return sorted([(graph.vertices[u[0]][u[1]], graph.vertices[v[0]][v[1]]) for u, v in graph.graph.edges])
#
#         return extract_edges(self) == extract_edges(other)
#
#     def create_subgraph(self, starting_vertices):
#         sub_nodes = set()
#         for layer_idx, vertex_idx in starting_vertices:
#             if layer_idx < len(self.vertices) and vertex_idx < len(self.vertices[layer_idx]):
#                 sub_nodes.add((layer_idx, vertex_idx))
#                 sub_nodes.update(nx.descendants(self.graph, (layer_idx, vertex_idx)))
#
#         subgraph_vertices = [[] for _ in range(len(self.vertices))]
#         for (layer_idx, vertex_idx) in sorted(sub_nodes):
#             if vertex_idx < len(self.vertices[layer_idx]):
#                 subgraph_vertices[layer_idx].append(self.vertices[layer_idx][vertex_idx])
#
#         subgraph_edges = [[] for _ in range(len(self.edges))]
#         for u, v in self.graph.edges:
#             if u in sub_nodes and v in sub_nodes:
#                 subgraph_edges[u[0]].append((subgraph_vertices[u[0]].index(self.vertices[u[0]][u[1]]),
#                                              subgraph_vertices[v[0]].index(self.vertices[v[0]][v[1]])))
#
#         return LayeredGraph(subgraph_vertices, subgraph_edges)
#
#     def create_quotient_graph(self, target_vertices):
#         quotient_nodes = set()
#         for layer_idx, vertex_idx in target_vertices:
#             if layer_idx < len(self.vertices) and vertex_idx < len(self.vertices[layer_idx]):
#                 quotient_nodes.add((layer_idx, vertex_idx))
#                 quotient_nodes.update(nx.ancestors(self.graph, (layer_idx, vertex_idx)))
#
#         quotient_vertices = [[] for _ in range(len(self.vertices))]
#         for (layer_idx, vertex_idx) in sorted(quotient_nodes):
#             if vertex_idx < len(self.vertices[layer_idx]):
#                 quotient_vertices[layer_idx].append(self.vertices[layer_idx][vertex_idx])
#
#         quotient_edges = [[] for _ in range(len(self.edges))]
#         for u, v in self.graph.edges:
#             if u in quotient_nodes and v in quotient_nodes:
#                 quotient_edges[u[0]].append((quotient_vertices[u[0]].index(self.vertices[u[0]][u[1]]),
#                                              quotient_vertices[v[0]].index(self.vertices[v[0]][v[1]])))
#
#         return LayeredGraph(quotient_vertices, quotient_edges)
#
#     def generate_all_subgraphs(self):
#         all_subgraphs = []
#         all_vertices = [(layer_idx, v_idx) for layer_idx, layer in enumerate(self.vertices) for v_idx in
#                         range(len(layer))]
#
#         for subset in range(1, len(all_vertices) + 1):
#             for sub_nodes in combinations(all_vertices, subset):
#                 subgraph = self.create_subgraph(sub_nodes)
#                 if not any(subgraph.is_isomorphic(existing) for existing in all_subgraphs):
#                     all_subgraphs.append(subgraph)
#         return all_subgraphs
#
#     def generate_all_quotient_graphs(self):
#         all_quotient_graphs = []
#         all_vertices = [(layer_idx, v_idx) for layer_idx, layer in enumerate(self.vertices) for v_idx in
#                         range(len(layer))]
#
#         for subset in range(1, len(all_vertices) + 1):
#             for target_nodes in combinations(all_vertices, subset):
#                 quotient_graph = self.create_quotient_graph(target_nodes)
#                 if not any(quotient_graph.is_isomorphic(existing) for existing in all_quotient_graphs):
#                     all_quotient_graphs.append(quotient_graph)
#         return all_quotient_graphs
#
#     def generate_functions(self, other):
#         functions = []
#         all_quotients = self.generate_all_quotient_graphs()
#         all_subgraphs = other.generate_all_subgraphs()
#
#         for q in all_quotients:
#             for s in all_subgraphs:
#                 if q.is_isomorphic(s):
#                     functions.append((q, s))
#         return functions
#
#     def generate_meta_graph(self):
#         subgraphs = self.generate_all_subgraphs()
#         meta_vertices = subgraphs
#         meta_edges = []
#
#         for i, sub1 in enumerate(subgraphs):
#             for j, sub2 in enumerate(subgraphs):
#                 if i != j and sub1.generate_functions(sub2):
#                     meta_edges.append((i, j))
#
#         return LayeredGraph(meta_vertices, meta_edges)
#
# # Example usage
# # Example usage
# if __name__ == "__main__":
#     graph3 = LayeredGraph([
#         ['A', 'B', 'C', 'A'],
#         ['B', 'C', 'A'],
#         ['A', 'B'],
#         ['C']
#     ], [
#         [(0, 0), (1, 0)],
#         [(1, 0), (2, 0), (1, 1)],
#         [(2, 0), (3, 0)]
#     ])
#
#     meta_graph = graph3.generate_meta_graph()
#     meta_graph.draw()
