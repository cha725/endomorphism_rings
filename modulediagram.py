import networkx as nx
from itertools import combinations, chain
import matplotlib.pyplot as plt
from typing import Optional
from submodule import Submodule
from quotientmodule import QuotientModule
import mod_graph_ops as mgo
from collections import defaultdict
from types import MappingProxyType

def all_sublists(n: int):
    L = list(range(n))
    S = [subset for subset in chain.from_iterable(combinations(L, r) for r in range(1,len(L) + 1))]
    S.sort(key=len)
    return S

def composition(m1 : dict, m2 : dict):
    return {k: m2[v] for k, v in m1.items() if v in m2}

class Arrow:
    def __init__(self, source : int, target : int, label : str):
        self.source = source
        self.target = target
        self.label = label

class ModuleDiagram:
    def __init__(self, 
                 arrows: list[Arrow] = [],
                 vertex_labels: Optional[tuple[str, ...]] = None):
        self.arrows = arrows

        G = nx.MultiDiGraph()
        if vertex_labels is not None:
            # Add vertices with labels
            for idx, vertex in enumerate(vertex_labels):
                G.add_node(idx, label=vertex)
        # Add edges with labels
        for arrow in self.arrows:
            G.add_edge(arrow.source, arrow.target, label=arrow.label)
        self.basic_graph = G
        if not nx.is_directed_acyclic_graph(self.basic_graph):
            raise ValueError("Invalid module diagram. Cannot contain a cycle.")
        self.vertex_labels = self.basic_graph.nodes.keys()
        self.num_vertices = len(self.vertex_labels)
        self.num_arrows = len(arrows)
        self._add_radical_labels()

    def _add_radical_labels(self):
        """
        Add a 'radical_layer' label to each node in the graph.
        Nodes in sources get layer 0, and each target node gets a layer
        one higher than the maximum of its predecessors.
        """
        # Process nodes in topological order
        sorted_nodes = nx.topological_sort(self.basic_graph)
        for node in sorted_nodes:
            # Layer is max of (predecessor layer + 1) or 0 if no predecessors
            pred_layer = []
            for pred in self.basic_graph.predecessors(node):
                pred_layer.append(self.basic_graph.nodes[pred].get("radical_layer",0))
            node_layer = max(pred_layer, default=-1) + 1
            self.basic_graph.nodes[node]["radical_layer"] = node_layer
    
    def nodes(self):
        return list(self.basic_graph.nodes)

    def node_to_radical_layers(self):
        """
        Returns:
            dict[int, int]: A dictionary from the node to its radical layer.
        """
        return {node : self.basic_graph.nodes[node]["radical_layer"] for node in self.basic_graph.nodes}

    def radical_layers_to_nodes(self):
        """
        Returns:
            dict[int, list[int]]: A dictionary from the radical layer to a list of nodes at that layer.
        """
        node_to_radical_layers = self.node_to_radical_layers()
        radical_layers_to_nodes = defaultdict(list)
        for node, layer in node_to_radical_layers.items():
            radical_layers_to_nodes[layer].append(node)
        return dict(radical_layers_to_nodes)
    
    def num_radical_layers(self):
        return len(self.radical_layers_to_nodes().keys())

    def draw_radical_layers(self):
        """
        Draw the nxgraph lined up with radical layers.
        """
        pos = {}
        for layer, nodes in self.radical_layers_to_nodes().items():
            n = len(nodes)
            if n == 1:
                x_coord = [0.5]
            else:
                x_coord = [i / (n - 1) for i in range(n)]
            for x, node in zip(x_coord, nodes):
                pos[node] = (x, -layer)
        nx.draw(self.basic_graph, pos, node_size=1000, node_color='lightblue', font_size=10,
                font_weight='bold')
        plt.show()

    def __eq__(self, other : ModuleDiagram):
        return list(nx.vf2pp_all_isomorphisms(self.basic_graph, other.basic_graph, node_label='label'))

    def node_bitmask(self):
        return MappingProxyType({node : 1 << idx for idx, node in enumerate(self.nodes())})

    def descendants_bitmask(self):
        """
        Need to be able to generate all submodules which is going to take forever.
        Idea: Create a bitmask that stores in position n answer to the question:
            Is descendant of node n?
        This way only have to search through the graph once and submodule creation
        comes down to bitwise operations.
        Especially since more than one generating set can create the same submodule.
        """
        bitmask = self.node_bitmask()
        descendants = {node : 0 for node in self.nodes()}
        for node in self.nodes():
            for des in nx.descendants(self.basic_graph, node):
                descendants[node] += bitmask[des]
        return MappingProxyType(descendants)
    
    def create_submodule_bitmask(self, gen_nodes : list[int]):
        """
        With the bitmasks should just take all descendant bitmasks and OR them.
        Any element that is a descendant of some generator is in the submodule.
        """
        node_bitmask = self.node_bitmask()
        des_bitmask = self.descendants_bitmask()
        sub_bitmask = 0
        for node in gen_nodes:
            sub_bitmask |= node_bitmask[node]
            sub_bitmask |= des_bitmask[node]
        return [node for node in self.nodes() if (node_bitmask[node] & sub_bitmask) == node_bitmask[node]]
    
    def generate_all_submodules_bitmask(self):
        node_bitmask = self.node_bitmask()
        des_bitmask = self.descendants_bitmask()
        submods = []
        for n in range(sum(node_bitmask.values())+1):
            submod_elts = [node for node, bitmask in node_bitmask.items() if (n & bitmask) == bitmask]
            if all((n & des_bitmask[node]) == des_bitmask[node] for node in submod_elts):
                submods.append([node for node, bitmask in node_bitmask.items() if (n & bitmask) == bitmask])
        return submods
            
    def ancestors_bitmask(self):
        """
        Same idea as submodule except now quotient module.
        """
        bitmask = self.node_bitmask()
        ancestors = {node : 0 for node in self.nodes()}
        for node in self.nodes():
            for des in nx.ancestors(self.basic_graph, node):
                ancestors[node] += bitmask[des]
        return MappingProxyType(ancestors)
    
    def create_quotientmodule_bitmask(self, gen_nodes : list[int]):
        """
        """
        node_bitmask = self.node_bitmask()
        anc_bitmask = self.ancestors_bitmask()
        sub_bitmask = 0
        for node in gen_nodes:
            sub_bitmask |= node_bitmask[node]
            sub_bitmask |= anc_bitmask[node]
        return [node for node in self.nodes() if (node_bitmask[node] & sub_bitmask) == node_bitmask[node]]
    
    def generate_all_quotientmodules_bitmask(self):
        node_bitmask = self.node_bitmask()
        anc_bitmask = self.ancestors_bitmask()
        submods = []
        for n in range(sum(node_bitmask.values())+1):
            submod_elts = [node for node, bitmask in node_bitmask.items() if (n & bitmask) == bitmask]
            if all((n & anc_bitmask[node]) == anc_bitmask[node] for node in submod_elts):
                submods.append([node for node, bitmask in node_bitmask.items() if (n & bitmask) == bitmask])
        return submods




    # def create_quotient_module(self, gen_nodes : list[int]):
    #     radical_layers_to_nodes = self.radical_layers_to_nodes()
    #     for layer in radical_layers_to_nodes.keys():
    #         layer_gen_nodes = [node for node in radical_layers_to_nodes[layer] if node in gen_nodes]
    #         predecessors = []
    #         for node in layer_gen_nodes:
    #             predecessors += list(self.basic_graph.predecessors(node))
    #         gen_nodes += predecessors
    #         gen_nodes = list(set(gen_nodes))
    #     return gen_nodes



    # def generate_all_quotient_modules(self):
    #     all_quotient_modules = [self.create_quotient_module([0])]
    #     all_generating_vertices = all_sublists(self.num_vertices)
    #     for generating_vertices in all_generating_vertices:
    #         new_quotient_module = self.create_quotient_module(generating_vertices)
    #         H = new_quotient_module[1].graph
    #         G = H.to_undirected()
    #         add = True
    #         while add:
    #             if nx.is_connected(G):
    #                 for quotient in all_quotient_modules:
    #                     if new_quotient_module[1].is_isomorphic_to(quotient[1]):
    #                         add = False
    #                 if add:
    #                     all_quotient_modules += [new_quotient_module]
    #                     add = False
    #             else:
    #                 add = False
    #     return all_quotient_modules

    # def homomorphism_group(self,other):
    #     hom = []
    #     for quotient in self.generate_all_quotient_modules():
    #         for sub in other.generate_all_sub_modules():
    #             for iso in quotient[1].is_isomorphic_to(sub[1]):
    #                 mapping = composition(quotient[2], composition(iso, sub[2]))
    #                 if mapping not in hom:
    #                     hom += [mapping]
    #     return hom


if __name__ == "__main__":


    arrows = [Arrow(1, 2, "b"), Arrow(2, 3, "c"), Arrow(4, 3, "d"), Arrow(0, 1, "a"), Arrow(0, 4, "e")]
    diagram = ModuleDiagram(arrows)

    print(diagram.basic_graph.nodes(data=True))
    print(diagram.basic_graph.edges(data=True))
    print(diagram.radical_layers_to_nodes())
    print("nodes", diagram.nodes())
    bit = diagram.node_bitmask()
    print("bitmask", bit)
    des = diagram.descendants_bitmask()
    print("descendants", des)



    print("\n All Submodules: ", diagram.generate_all_submodules_bitmask())

#     MB = modulediagram(['2', '2', '1', '3', '1', '3'], [(1, 2), (1, 3), (1, 4), (2, 5), (3, 5), (3, 6), (4, 6)],[1,2,3,4,5,6])

#     #print('vertices', MB.graph.nodes())
#     #print('radical layers', MB.radical_layers)
#     #print(MB.create_radical_labels())
#     print(MB.create_submodule([2,2,3]).sub_vert)
#     #print(MB.generate_all_sub_modules())
#     #print(MB.graph.nodes())
# #
# #     G = nx.DiGraph()
# #     G.add_edges_from([(0,1),(1,2),(2,3),(3,4),(4,5)])
# #     print(source_induced_subgraph(G,[1]))
# #     print(sink_induced_subgraph(G,[1]))


# #### Old code

#     #
#     #
#     # def create_radical_layers(self):
#     #     return mgo.create_radical_layers(self.basic_graph)
#     #
#     # # def create_radical_layers(self):  # G is a graph in networkx
#     # #     # A vertex is in radical layer l if it is l steps away from a vertex that is a source
#     # #     G = self.basic_graph
#     # #     sources = [x for x in G.nodes() if G.in_degree(x) == 0]
#     # #     if sources == []:
#     # #         print('This nx.DiGraph has no vertices that are sources.')
#     # #         return False
#     # #     radical_layers = [sources]
#     # #     test_vertices = sources
#     # #     vertices_left = list(G.nodes())
#     # #     while vertices_left != []:  # This is to make sure each vertex is assigned a unique radical level
#     # #         next_neighbors = []
#     # #         for vertex in test_vertices:
#     # #             next_neighbors += [x for x in G.neighbors(vertex) if x in vertices_left]
#     # #             if vertex in vertices_left:
#     # #                 vertices_left.remove(vertex)
#     # #         if next_neighbors != []:
#     # #             radical_layers += [next_neighbors]
#     # #             test_vertices = next_neighbors
#     # #         else:
#     # #             return radical_layers
#     #
#     # # def create_radical_labels(self):
#     # #     new_labels = []
#     # #     for layer in range(len(self.radical_layers)):
#     # #         vertices = self.radical_layers[layer]
#     # #         print('vert', vertices)
#     # #         print('mapping', self.indices_to_labels)
#     # #         new_labels.append([self.indices_to_labels[vertex] for vertex in vertices])
#     # #     return new_labels
#     #
#     # def find_top(self):
#     #     # The top of a DiGraph is the set of vertices that are sources
#     #     if self.radical_layers[0]:
#     #         return self.radical_layers[0]
#     #     else:
#     #         print('This nxDiGraph has no vertices that are sources.')
#     #         return False
#     #
#     # def find_socle(self):
#     #     # The socle of a DiGraph is the set of vertices that are sinks
#     #     G = self.basic_graph
#     #     socle = [x for x in G.nodes() if G.out_degree(x) == 0]
#     #     if socle == []:
#     #         print('This nx.DiGraph has no vertices that are sinks.')
#     #         return False
#     #     else:
#     #         return socle