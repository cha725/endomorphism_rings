import networkx as nx
from itertools import combinations, chain
import matplotlib.pyplot as plt

from history import ModuleDiagram
from submodule import Submodule
from quotientmodule import QuotientModule
import mod_graph_ops as mgo

def all_sublists(n: int):
    L = list(range(n))
    S = [subset for subset in chain.from_iterable(combinations(L, r) for r in range(1,len(L) + 1))]
    S.sort(key=len)
    return S

def composition(m1 : dict, m2 : dict):
    return {k: m2[v] for k, v in m1.items() if v in m2}


class modulediagram:
    def __init__(self, vertex_labels, arrows, vertex_indices=None, arrow_labels=None):
        self.vertex_labels = vertex_labels # list of strings; one for each vertex
        self.num_vertices = len(vertex_labels)
        self.arrows = arrows  # list of tuples of vertex indices; (source,sink)
        self.num_arrows = len(arrows)
        self.vertex_indices = vertex_indices
        self.arrow_labels = arrow_labels

        self.indices = self.new_indices()
        self.indices_to_labels = self.vertex_index_to_vertex_label()

        self.basic_graph = self._create_basic_graph()
        self.radical_layers = mgo.create_radical_layers(self.basic_graph)
        self.graph = mgo.create_layered_graph(self.basic_graph)
        self.top = mgo.sources(self.basic_graph)
        self.socle = mgo.sinks(self.basic_graph)

        #self.all_submodules = self.generate_all_sub_modules()
        #self.all_quotient_modules = self.generate_all_quotient_modules()

    def new_indices(self):
        if self.vertex_indices is None:
            new_indices = range(self.num_vertices)
        else:
            new_indices = self.vertex_indices
        return new_indices

    def vertex_index_to_vertex_label(self):
        mapping = {}
        for index_id, index in enumerate(self.new_indices()):
            mapping[index] = self.vertex_labels[index_id]
        return mapping

    def _create_basic_graph(self):
        g = nx.MultiDiGraph()
        for i in self.new_indices():
            g.add_node(i, label = self.indices_to_labels[i])
        g.add_edges_from(self.arrows)
        return g


    #
    # def relabel_mapping_wrt_radical_layers(self):
    #     # This changes the labelling of the vertices in the DiGraph so that the vertices in radical layer l are strictly smaller than those in radical layer l+1
    #     mapping_canonical_to_example = {}
    #     mapping_example_to_canonical = {}
    #     counter = 0
    #     for layer in self.radical_layers:
    #         for vertex in self.radical_layers[layer]:
    #             mapping_canonical_to_example[counter] = vertex
    #             mapping_example_to_canonical[vertex] = counter
    #             counter += 1
    #     return mapping_canonical_to_example, mapping_example_to_canonical

    def _create_layered_graph(self): # takes the ModuleDiagram and interprets it as a networkx graph
        G = nx.DiGraph()
        for layer in self.radical_layers:
            for vertex in self.radical_layers[layer]:
                G.add_node(vertex, layer=layer, label= self.indices_to_labels[vertex])
        for arrow in self.arrows:
            G.add_edge(*arrow)
        return G

    def draw(self): # draws the networkx graph nicely
        pos = {}
        labels = {}
        for layer in self.radical_layers:
            for vertex_idx, vertex in enumerate(self.radical_layers[layer]):
                pos[vertex] = (vertex_idx, -layer)
                labels[vertex] = vertex,self.indices_to_labels[vertex]

        nx.draw(self.graph, pos, labels=labels, with_labels=True, node_size=1000, node_color='lightblue', font_size=10,
                font_weight='bold', arrows=True)
        plt.show()

    def create_submodule(self, gens: list):
        g = mgo.out_directed_subgraph(self.graph, gens)
        print(g.nodes(data=True))
        inclusion = mgo.relabelling_map_wrt_radical_layers(g)
        print(inclusion)
        sub_g = mgo.relabelling_wrt_radical_layers(g)
        print(sub_g.nodes(data=True))
        return Submodule(self.vertex_labels, self.arrows, list(sub_g.nodes()), inclusion)

    def create_quotient_module(self, gens:list):
        g = mgo.in_directed_subgraph(self.graph, gens)
        projection = mgo.relabelling_map_wrt_radical_layers(g)
        quot_g = mgo.relabelling_wrt_radical_layers(g)
        return QuotientModule(self.vertex_labels, self.arrows, list(quot_g.nodes()), projection)

    # def create_quotient_module(self, generating_vertices):
    #     if not generating_vertices:
    #         return
    #     generating_vertices = list(set(generating_vertices))
    #     test_vertices = generating_vertices
    #     all_vertices = generating_vertices
    #     tested_vertices = []
    #     while test_vertices:
    #         new_vertices = []
    #         for vertex in test_vertices:
    #             if vertex not in tested_vertices:
    #                 new_vertices += list(self.graph.predecessors(vertex))
    #                 if vertex not in all_vertices:
    #                     all_vertices += [vertex]
    #                 tested_vertices += [vertex]
    #         test_vertices = new_vertices
    #     H = nx.induced_subgraph(self.graph, all_vertices)
    #     mapping = {}
    #     for vertex_id, vertex in enumerate(H.nodes()):
    #         mapping[vertex] = vertex_id
    #     HH = nx.relabel_nodes(H, mapping)
    #     vertex_labels = [HH.nodes(data=True)[i]['label'] for i in HH.nodes()]
    #     arrows = [a for a in HH.edges]
    #     M = ModuleDiagram(vertex_labels, arrows)
    #     return [self, M, mapping]


    def is_isomorphic_to(self,other):
        # def node_match(n1, n2):
        #     return n1['layer'] == n2['layer'] and n1['label'] == n2['label']
        # nx.is_isomorphic(self.graph, other.graph, node_match = node_match),
        return list(nx.vf2pp_all_isomorphisms(self.graph, other.graph, node_label='label'))

    def generate_all_sub_modules(self):
        all_submodules = [self.create_submodule([0])]
        all_generating_vertices = all_sublists(self.num_vertices)
        for generating_vertices in all_generating_vertices:
            new_submodule = self.create_submodule(generating_vertices)
            submodule_diagram = ModuleDiagram()
            H = new_submodule[1].graph
            G = H.to_undirected()
            add = True
            while add:
                if nx.is_connected(G):
                    for submodule in all_submodules:
                        if new_submodule[1].is_isomorphic_to(submodule[1]):
                            add = False
                            break
                    if add:
                        all_submodules += [new_submodule]
                        add = False
                else:
                    add = False
        return all_submodules

    def generate_all_quotient_modules(self):
        all_quotient_modules = [self.create_quotient_module([0])]
        all_generating_vertices = all_sublists(self.num_vertices)
        for generating_vertices in all_generating_vertices:
            new_quotient_module = self.create_quotient_module(generating_vertices)
            H = new_quotient_module[1].graph
            G = H.to_undirected()
            add = True
            while add:
                if nx.is_connected(G):
                    for quotient in all_quotient_modules:
                        if new_quotient_module[1].is_isomorphic_to(quotient[1]):
                            add = False
                    if add:
                        all_quotient_modules += [new_quotient_module]
                        add = False
                else:
                    add = False
        return all_quotient_modules

    def homomorphism_group(self,other):
        hom = []
        for quotient in self.generate_all_quotient_modules():
            for sub in other.generate_all_sub_modules():
                for iso in quotient[1].is_isomorphic_to(sub[1]):
                    mapping = composition(quotient[2], composition(iso, sub[2]))
                    if mapping not in hom:
                        hom += [mapping]
        return hom


if __name__ == "__main__":

    MB = modulediagram(['2', '2', '1', '3', '1', '3'], [(1, 2), (1, 3), (1, 4), (2, 5), (3, 5), (3, 6), (4, 6)],[1,2,3,4,5,6])

    #print('vertices', MB.graph.nodes())
    #print('radical layers', MB.radical_layers)
    #print(MB.create_radical_labels())
    print(MB.create_submodule([2,2,3]).sub_vert)
    #print(MB.generate_all_sub_modules())
    #print(MB.graph.nodes())
#
#     G = nx.DiGraph()
#     G.add_edges_from([(0,1),(1,2),(2,3),(3,4),(4,5)])
#     print(source_induced_subgraph(G,[1]))
#     print(sink_induced_subgraph(G,[1]))


#### Old code

    #
    #
    # def create_radical_layers(self):
    #     return mgo.create_radical_layers(self.basic_graph)
    #
    # # def create_radical_layers(self):  # G is a graph in networkx
    # #     # A vertex is in radical layer l if it is l steps away from a vertex that is a source
    # #     G = self.basic_graph
    # #     sources = [x for x in G.nodes() if G.in_degree(x) == 0]
    # #     if sources == []:
    # #         print('This nx.DiGraph has no vertices that are sources.')
    # #         return False
    # #     radical_layers = [sources]
    # #     test_vertices = sources
    # #     vertices_left = list(G.nodes())
    # #     while vertices_left != []:  # This is to make sure each vertex is assigned a unique radical level
    # #         next_neighbors = []
    # #         for vertex in test_vertices:
    # #             next_neighbors += [x for x in G.neighbors(vertex) if x in vertices_left]
    # #             if vertex in vertices_left:
    # #                 vertices_left.remove(vertex)
    # #         if next_neighbors != []:
    # #             radical_layers += [next_neighbors]
    # #             test_vertices = next_neighbors
    # #         else:
    # #             return radical_layers
    #
    # # def create_radical_labels(self):
    # #     new_labels = []
    # #     for layer in range(len(self.radical_layers)):
    # #         vertices = self.radical_layers[layer]
    # #         print('vert', vertices)
    # #         print('mapping', self.indices_to_labels)
    # #         new_labels.append([self.indices_to_labels[vertex] for vertex in vertices])
    # #     return new_labels
    #
    # def find_top(self):
    #     # The top of a DiGraph is the set of vertices that are sources
    #     if self.radical_layers[0]:
    #         return self.radical_layers[0]
    #     else:
    #         print('This nxDiGraph has no vertices that are sources.')
    #         return False
    #
    # def find_socle(self):
    #     # The socle of a DiGraph is the set of vertices that are sinks
    #     G = self.basic_graph
    #     socle = [x for x in G.nodes() if G.out_degree(x) == 0]
    #     if socle == []:
    #         print('This nx.DiGraph has no vertices that are sinks.')
    #         return False
    #     else:
    #         return socle