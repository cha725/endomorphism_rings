import networkx as nx
import matplotlib.pyplot as plt
from typing import Optional
from collections import defaultdict
from functools import cached_property
from quiver_algebra import Arrow
from bitmask_subgraph import BitmaskSubgraph


class ModuleDiagram:
    def __init__(self, 
                 arrows: list[Arrow] = [],
                 isolated_vertices : Optional[list[int]] = None,
                 vertex_simples : Optional[dict[int,int]] = None,
                 vertex_labels: Optional[dict[int,str]] = None):
        self.arrows = arrows
        self.bitmaskgraph = BitmaskGraph(self.arrows)
        self.vertex_labels = vertex_labels

        G = nx.MultiDiGraph()
        if self.vertex_labels is not None:
            # Add vertices with labels
            for vertex, label in self.vertex_labels.items():
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
        self.vertex_simples = vertex_simples
        self.num_vertices = len(self.basic_graph.nodes)
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

        nx.draw(self.basic_graph, pos,
            node_size=1000, node_color='lightblue',
            font_weight='bold')
        nx.draw_networkx_labels(self.basic_graph, pos, labels=self.vertex_simples,
                            font_size=10, font_weight='bold')  
        edge_labels = {
            (u, v, k): d.get("label")
            for u, v, k, d in self.basic_graph.edges(keys=True, data=True)
            if d.get("label") is not None
        }

        nx.draw_networkx_edge_labels(
            self.basic_graph, pos,
            edge_labels=edge_labels,
            font_size=9,
            label_pos=0.5,
            rotate=False
        )      
        plt.show()

    def __eq__(self, other : ModuleDiagram):
        def node_match(node1_att,node2_att):
            return node1_att.get("simples") == node2_att.get("simples")
        def edge_match(edge1_att,edge2_att):
            return edge1_att.get("label") == edge2_att.get("label")
        return nx.is_isomorphic(self.basic_graph,other.basic_graph,node_match=node_match,edge_match=edge_match)
    
    @cached_property
    def generate_all_submodules(self):
        return self.bitmaskgraph.compute_succ_closed_subsets
    
    @cached_property
    def generate_all_quotients(self):
        return self.bitmaskgraph.compute_pred_closed_subsets
    
    def is_submodule(self, other : ModuleDiagram):
        return other in self.generate_all_submodules
    
    def is_quotient(self, other : ModuleDiagram):
        return other in self.generate_all_quotients
        
    def hom_group(self, other : ModuleDiagram):
        hom = []
        for quotient in self.generate_all_quotients:
            quotdiagram = QuotientModuleDiagram(self,quotient)
            for submod in other.generate_all_submodules:
                submoddiagram = SubModuleDiagram(other,submod)
                if quotdiagram == submoddiagram:
                    hom.append((quotdiagram, submoddiagram))
        return hom
    
    def __repr__(self):
        if not self.nodes:
            return f"Module diagram with no vertices and no arrows."
        if not self.arrows:
            return f"Module diagram with vert={self.vertex_labels} and no arrows."
        return f"Module diagram with vert={self.vertex_labels} and arrows = {self.arrows}"

class QuotientModuleDiagram(ModuleDiagram):
    def __init__(self,
                 parent : ModuleDiagram,
                 vertices : list[int]):
        if not set(vertices) <= set(parent.nodes):
            raise ValueError(f"Invalid vertices. {vertices} must be a subset of {parent.nodes}")
        self.parent = parent
        subgraph = parent.basic_graph.subgraph(vertices).copy()
        arrows = [Arrow(u, v, d.get("label")) for u, v, d in subgraph.edges(data=True)]
        vertex_labels : dict = {n: subgraph.nodes[n].get("label") for n in subgraph.nodes}
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
        vertex_labels :dict = {n: subgraph.nodes[n].get("label") for n in subgraph.nodes}
        vertex_simples = {n: parent.basic_graph.nodes[n]["simples"] for n in subgraph.nodes}
        super().__init__(arrows=arrows, vertex_labels=vertex_labels, vertex_simples=vertex_simples)




class Examples:
    """
    Class to store examples of module diagrams.
    """
    def __init__(self, 
                 examples: dict[str, dict]):
        self.examples = examples

    def add(self, 
            name: str, 
            arrows: list, 
            isolated_vertices: Optional[list] = None, 
            vertex_simples: Optional[dict] = None):
        self.examples[name] = {
            "arrows": arrows,
            "isolated_vertices": isolated_vertices or [],
            "vertex_simples": vertex_simples or {}
        }

    def run(self, draw: bool = False):
        for name, data in self.examples.items():
            print(f"\n=== Example: {name} ===")
            diagram = ModuleDiagram(arrows=data["arrows"],
                                    isolated_vertices=data["isolated_vertices"],
                                    vertex_simples=data["vertex_simples"])

            if draw:
                diagram.draw_radical_layers

            print("\nNodes:", diagram.nodes)
            print("\nRadical layers:", diagram.node_to_radical_layers)
            print("\nAll submodules:", diagram.generate_all_submodules)
            print("\nAll quotient modules:", diagram.generate_all_quotients)
            # print("\nEndomorphisms:")
            # for idx, hom in enumerate(diagram.hom_group(diagram)):
            #     print(f"{idx}: {hom}")

# TODO: make the print out prettier.

if __name__ == "__main__":

    examples = Examples({})

    examples.add("Example 1: A4",
                        arrows=[Arrow(0,1,"a"), Arrow(1,2,"b"), Arrow(2,3,"c")],
                        vertex_simples={0:0, 1:1, 2:0, 3:1})

    examples.add("Example 2: Diamond",
                        arrows=[Arrow(0,1,"a"), Arrow(0,2,"b"), Arrow(1,3,"c"), Arrow(2,3,"d")])

    examples.add("Example 3: Y structure",
                        arrows=[Arrow(0,2,"a"), Arrow(1,2,"b"), Arrow(2,3,"c"), Arrow(3,4,"d")])

    examples.add("Example 4: weird radical",
                        arrows=[Arrow(0,1,"a"), Arrow(1,2,"b"), Arrow(2,3,"c"), Arrow(3,4,"d"),
                                Arrow(0,5,"e"), Arrow(5,4,"f")])

    # examples.add("Example 5: product of modules",
    #                     arrows=[Arrow(0,1,"a1"), Arrow(1,2,"b1"), Arrow(2,3,"c1"),
    #                             Arrow(4,5,"a2"), Arrow(5,6,"b2"), Arrow(6,7,"c2")])

    # examples.add("Example 6: Algebra and dual",
    #                     arrows=[Arrow(0,1,"a1"), Arrow(0,2,"b1"), Arrow(4,5,"a2"), Arrow(6,7,"b2")],
    #                     isolated_vertices=[3])

    # examples.add("Example 7: Random graph",
    #                     arrows=[],  
    #                     isolated_vertices=[])
    
    examples.run()

