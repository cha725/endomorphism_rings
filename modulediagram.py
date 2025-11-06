import networkx as nx
import matplotlib.pyplot as plt
from typing import Optional
from collections import defaultdict
from functools import cached_property
from quiver_algebra import Arrow
from bitmaskgraph import BitmaskGraph


class ModuleDiagram:
    def __init__(self, 
                 arrows: list[Arrow] = [],
                 isolated_vertices : Optional[list[int]] = None,
                 vertex_simples : Optional[dict[int,int]] = None,
                 vertex_labels: Optional[dict[int,str]] = None):
        self.arrows = arrows
        self.bitmaskgraph = BitmaskGraph(self.arrows)

        G = nx.MultiDiGraph()
        if vertex_labels is not None:
            # Add vertices with labels
            for vertex, label in vertex_labels.items():
                G.add_node(vertex, label=label)

        for arrow in self.arrows:
            G.add_edge(arrow.source, arrow.target, label=arrow.label)

        if isolated_vertices is not None:
            for vertex in isolated_vertices:
                G.add_node(vertex)
        self.basic_graph = G
        if not nx.is_directed_acyclic_graph(self.basic_graph):
            raise ValueError("Invalid module diagram. Cannot contain a cycle.")
        for vertex in self.basic_graph.nodes:
            self.basic_graph.nodes[vertex]["simples"] = vertex_simples.get(vertex) if vertex_simples else None
        self.vertex_labels = self.basic_graph.nodes.keys()
        self.num_vertices = len(self.vertex_labels)
        self.num_arrows = len(arrows)
        self._add_radical_labels()
        self.nodes = list(nx.topological_sort(self.basic_graph))

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
    
    @cached_property
    def node_to_radical_layers(self):
        """
        Returns:
            dict[int, int]: A dictionary from the node to its radical layer.
        """
        return {node : self.basic_graph.nodes[node]["radical_layer"] for node in self.basic_graph.nodes}

    @cached_property
    def radical_layers_to_nodes(self):
        """
        Returns:
            dict[int, list[int]]: A dictionary from the radical layer to a list of nodes at that layer.
        """
        node_to_radical_layers = self.node_to_radical_layers
        radical_layers_to_nodes = defaultdict(list)
        for node, layer in node_to_radical_layers.items():
            radical_layers_to_nodes[layer].append(node)
        return dict(radical_layers_to_nodes)
    
    def ordered_nodes_wrt_radical(self):
        radical_layers_to_nodes = self.radical_layers_to_nodes
        ordered = []
        for layer in range(len(radical_layers_to_nodes.items())):
            ordered += radical_layers_to_nodes[layer]
        return ordered
    
    def num_radical_layers(self):
        return len(self.radical_layers_to_nodes.keys())

    @cached_property
    def draw_radical_layers(self):
        """
        Draw the nxgraph lined up with radical layers.
        """
        pos = {}
        for layer, nodes in self.radical_layers_to_nodes.items():
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
        def node_match(node1_att,node2_att):
            return node1_att.get("simples") == node2_att.get("simples")
        def edge_match(edge1_att,edge2_att):
            return edge1_att.get("label") == edge2_att.get("label")
        return nx.is_isomorphic(self.basic_graph,other.basic_graph,node_match=node_match,edge_match=edge_match)
    

    def generate_all_submodules(self):
        return self.bitmaskgraph.des_closed
    
    def generate_all_quotients(self):
        return self.bitmaskgraph.anc_closed
    
    # def generate_all_quotients_from_submod(self):
    #     """
    #     speeds up by an average of 35.49=time1/time2
    #     compared to generate all quotients
    #     """
    #     return self.bitmaskgraph.anc_closed_from_des
        


    
    # def create_submodule_bitmask(self, gen_nodes : list[int]):
    #     """
    #     With the bitmasks should just take all descendant bitmasks and OR them.
    #     Any element that is a descendant of some generator is in the submodule.
    #     """
    #     node_bitmask = self.node_bitmask
    #     des_bitmask = self.descendants_bitmask
    #     sub_bitmask = 0
    #     for node in gen_nodes:
    #         sub_bitmask |= node_bitmask[node]
    #         sub_bitmask |= des_bitmask[node]
    #     return [node for node in self.nodes if (node_bitmask[node] & sub_bitmask) == node_bitmask[node]]
    
    # def is_submodule(self, elements : list[int]):
    #     node_bitmask = self.node_bitmask
    #     des_bitmask = self.descendants_bitmask
    #     submod_bitmask = sum([node_bitmask[element] for element in elements])
    #     return all((submod_bitmask & des_bitmask[node]) == des_bitmask[node] for node in elements)

    # @cached_property
    # def generate_all_submodules_bitmask(self):
    #     """
    #     This generates all submodules not just indecomposable ones.
    #     """
    #     node_bitmask = self.node_bitmask
    #     total_node_bitmask = sum(node_bitmask.values())
    #     des_bitmask = self.descendants_bitmask
    #     submods = []
    #     for n in range(total_node_bitmask+1):
            


    #         submod_elts = [node for node, bitmask in node_bitmask.items() if (n & bitmask) == bitmask]
    #         if all((n & des_bitmask[node]) == des_bitmask[node] for node in submod_elts):
    #             submods.append(n)
    #     return submods
        
    # def bitmask_to_nodes(self, bitmask_list: list[int]):
    #     """
    #     Converts a list of submodule bitmasks to lists of nodes.
    #     """
    #     result = []
    #     for bitmask in bitmask_list:
    #         nodes = [node for node, value in self.node_bitmask.items() if (bitmask & value) == value]
    #         result.append(nodes)
    #     return result

    # def generate_all_submodules(self):
    #     return self.bitmask_to_nodes(self.generate_all_submodules_bitmask)

    # def create_quotientmodule_bitmask(self, gen_nodes : list[int]):
    #     """
    #     For every subgraph closed under predecessors, the complementary subgraph
    #     is closed under successors.
    #     So finding all quotient modules comes down to taking the complementary
    #     set for each subgraph.
    #     Since generate all subgraphs is a cached property doing this will not recompute
    #     the subgraphs.
    #     """
    #     node_bitmask = self.node_bitmask
    #     anc_bitmask = self.ancestors_bitmask
    #     sub_bitmask = 0
    #     for node in gen_nodes:
    #         sub_bitmask |= node_bitmask[node]
    #         sub_bitmask |= anc_bitmask[node]
    #     return [node for node in self.nodes if (node_bitmask[node] & sub_bitmask) == node_bitmask[node]]
    
    # @cached_property
    # def generate_all_quotient_modules_bitmask(self):
    #     quotmods = []
    #     n = self.num_vertices
    #     for submod in self.generate_all_submodules_bitmask:
    #         quotmods.append(submod ^ (2**n-1))
    #     return quotmods
    
    # def generate_all_quotient_modules(self):
    #     return self.bitmask_to_nodes(self.generate_all_quotient_modules_bitmask)
    
    def hom_group(self, other: ModuleDiagram):
        hom = []
        for quotient in self.generate_all_quotients():
            quotdiagram = QuotientModuleDiagram(self,quotient)
            for submod in other.generate_all_submodules():
                submoddiagram = SubModuleDiagram(other,submod)
                if quotdiagram == submoddiagram:
                    hom.append((quotdiagram, submoddiagram))
        return hom
    
    def __repr__(self):
        if not self.nodes:
            return f"Module diagram with no vertices and no arrows."
        if not self.arrows:
            return f"Module diagram with vertices {self.nodes} and no arrows."
        return f"Module diagram with vertices {self.nodes} and arrows {self.arrows}."
    
    

class QuotientModuleDiagram(ModuleDiagram):
    def __init__(self,
                 parent : ModuleDiagram,
                 vertices : list[int]):
        if not set(vertices) <= set(parent.nodes):
            raise ValueError(f"Invalid vertices. {vertices} must be a subset of {parent.nodes}")
        self.parent = parent
        subgraph = parent.basic_graph.subgraph(vertices).copy()
        arrows = [Arrow(u, v, d.get("label")) for u, v, d in subgraph.edges(data=True)]
        vertex_labels = {n: subgraph.nodes[n].get("label") for n in subgraph.nodes}
        vertex_simples = {n: parent.basic_graph.nodes[n]["simples"] for n in subgraph.nodes}
        super().__init__(arrows=arrows, vertex_labels=vertex_labels, vertex_simples=vertex_simples)


class SubModuleDiagram(ModuleDiagram):
    def __init__(self,
                 parent : ModuleDiagram,
                 vertices : list[int]):
        if not set(vertices) <= set(parent.nodes):
            raise ValueError(f"Invalid vertices. Must be a subset of {parent.nodes}")
        self.parent = parent
        subgraph = parent.basic_graph.subgraph(vertices).copy()
        arrows = [Arrow(u, v, d.get("label")) for u, v, d in subgraph.edges(data=True)]
        vertex_labels = {n: subgraph.nodes[n].get("label") for n in subgraph.nodes}
        vertex_simples = {n: parent.basic_graph.nodes[n]["simples"] for n in subgraph.nodes}
        super().__init__(arrows=arrows, vertex_labels=vertex_labels, vertex_simples=vertex_simples)
    


class ProjectiveModuleDiagram(ModuleDiagram):
    def __init__(self,
                 quiver_algebra,
                 vertex):
        pass



if __name__ == "__main__":
    import random, time
    examples = [1,2,3,4,5,6]

    examples = {  
        "Example 1: A4": {  
            "enabled": 1 in examples,  
            "arrows": [Arrow(0, 1, "a"), Arrow(1, 2, "b"), Arrow(2, 3, "c")],  
            "isolated_vertices": [],
            "simples": {0:0,1:1,2:0,3:1}
        },  
        "Example 2: k[x,y] diamond": {  
            "enabled": 2 in examples,  
            "arrows": [Arrow(0, 1, "a"), Arrow(0, 2, "b"), Arrow(1, 3, "c"), Arrow(2, 3, "d")],  
            "isolated_vertices": []  
        },  
        "Example 3: Y structure": {  
            "enabled": 3 in examples,  
            "arrows": [Arrow(0, 2, "a"), Arrow(1, 2, "b"), Arrow(2, 3, "c"), Arrow(3, 4, "d")],  
            "isolated_vertices": []  
        },  
        "Example 4: weird radical": {  
            "enabled": 4 in examples,  
            "arrows": [Arrow(0, 1, "a"), Arrow(1, 2, "b"), Arrow(2, 3, "c"), Arrow(3, 4, "d"), Arrow(0, 5, "e"), Arrow(5, 4, "f")],  
            "isolated_vertices": []  
        },  
        "Example 5: product of modules": {  
            "enabled": 5 in examples,  
            "arrows": [Arrow(0, 1, "a1"), Arrow(1, 2, "b1"), Arrow(2, 3, "c1"), Arrow(4, 5, "a2"), Arrow(5, 6, "b2"), Arrow(6, 7, "c2")],  
            "isolated_vertices": []  
        },  
        "Example 6: Algebra and dual": {  
            "enabled": 6 in examples,  
            "arrows": [Arrow(0, 1, "a1"), Arrow(0, 2, "b1"), Arrow(4, 5, "a2"), Arrow(6, 7, "b2")],  
            "isolated_vertices": [3]  
        },  
        "Example 7: Random graph": {  
            "enabled": 7 in examples,  
            "arrows": [],  # Generated below  
            "isolated_vertices": []  
        }  
    }  

    draw = False

  
    time_ratio = []
    rounds = 1
    for round in range(rounds):
        if examples["Example 7: Random graph"]["enabled"]:  
            n = 5
            for i in range(n):  
                for j in range(i+1, n):  
                    if random.random() < 0.2:  
                        examples["Example 7: Random graph"]["arrows"].append(Arrow(i, j))  

        for name, data in examples.items():  
            if not data["enabled"]:  
                continue  
            print(round)
            print(f"\n--- {name} ---")  
            diagram = ModuleDiagram(data["arrows"],
                                    isolated_vertices=data.get("isolated_vertices", []),
                                    vertex_simples=data.get("simples", []))  

            if draw:  
                diagram.draw_radical_layers  

            start = time.time()  
            print("\nNodes:", diagram.nodes)  
            print("\nRadical layers:", diagram.node_to_radical_layers)  
            print("\nAll submodules:", diagram.generate_all_submodules())  
            print("\nAll quotient modules:", diagram.generate_all_quotients())
            print("\nEndomorphism groups:")
            for idx, hom in enumerate(diagram.hom_group(diagram)):
                print(f"{idx}: {hom}")
            end = time.time()  
            print(f"\nCompleted in {end - start:.4f} seconds.")  

