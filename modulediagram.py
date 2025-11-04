import networkx as nx
import matplotlib.pyplot as plt
from typing import Optional
from collections import defaultdict
from types import MappingProxyType

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
    
    def ordered_nodes_wrt_radical(self):
        radical_layers_to_nodes = self.radical_layers_to_nodes()
        ordered = []
        for layer in range(len(radical_layers_to_nodes.items())):
            ordered += radical_layers_to_nodes[layer]
        return ordered
    
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
        bitmask = {}
        nodes = self.ordered_nodes_wrt_radical()
        for idx, node in enumerate(nodes):
            self.basic_graph.nodes[node]["bitmask"] = 1 << idx
            bitmask[node] = 1 << idx
        return MappingProxyType(bitmask)

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
        des_bitmask = {k : v for k, v in bitmask.items()}
        radical_layers_to_nodes = self.radical_layers_to_nodes()
        for _, node_list in sorted(radical_layers_to_nodes.items(), reverse=True):
            for node in node_list:
                for des in self.basic_graph.successors(node):
                    des_bitmask[node] |= des_bitmask[des]
        return MappingProxyType(des_bitmask)
    
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
    
    def is_submodule(self, elements : list[int]):
        node_bitmask = self.node_bitmask()
        des_bitmask = self.descendants_bitmask()
        submod_bitmask = sum([node_bitmask[element] for element in elements])
        return all((submod_bitmask & des_bitmask[node]) == des_bitmask[node] for node in elements)

    def generate_all_submodules_bitmask(self):
        """
        This generates all submodules not just indecomposable ones.
        """
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
        """
        bitmask = self.node_bitmask()
        anc_bitmask = {k : v for k, v in bitmask.items()}
        radical_layers_to_nodes = self.radical_layers_to_nodes()
        for _, node_list in sorted(radical_layers_to_nodes.items()):
            for node in node_list:
                for pred in self.basic_graph.predecessors(node):
                    anc_bitmask[node] |= anc_bitmask[pred]
        return MappingProxyType(anc_bitmask)
    
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
    


if __name__ == "__main__":

    # Example 1: A4 0 -> 1 -> 2 -> 3
    print("\n--- Example 1: A4 ---")
    arrows = [
        Arrow(0, 1, "a"),
        Arrow(1, 2, "b"),
        Arrow(2, 3, "c")
    ]
    diagram1 = ModuleDiagram(arrows)
    diagram1.draw_radical_layers()

    print("Nodes:", diagram1.nodes())
    print("Radical layers:", diagram1.node_to_radical_layers())
    print("All submodules:", diagram1.generate_all_submodules_bitmask())
    print("All quotient modules:", diagram1.generate_all_quotientmodules_bitmask())

    # Example 2: k[x,y] relations = {x^2=y^2, xy, yx}
    print("\n--- Example 2: k[x,y] diamond ---")
    arrows2 = [
        Arrow(0, 1, "a"),
        Arrow(0, 2, "b"),
        Arrow(1, 3, "c"),
        Arrow(2, 3, "d")
    ]
    diagram2 = ModuleDiagram(arrows2)
    diagram2.draw_radical_layers()

    print("Nodes:", diagram2.nodes())
    print("Radical layers:", diagram2.node_to_radical_layers())
    print("All submodules:", diagram2.generate_all_submodules_bitmask())
    print("All quotient modules:", diagram2.generate_all_quotientmodules_bitmask())


    # Example 3: Y structure
    print("\n--- Example 3: X structure ---")
    arrows3 = [
        Arrow(0, 2, "a"),
        Arrow(1, 2, "b"),
        Arrow(2, 3, "c"),
        Arrow(3, 4, "d"),
        Arrow(3, 5, "e")
    ]
    diagram3 = ModuleDiagram(arrows3)
    diagram3.draw_radical_layers()

    print("Nodes:", diagram3.nodes())
    print("Radical layers:", diagram3.node_to_radical_layers())
    print("All submodules:", diagram3.generate_all_submodules_bitmask())
    print("All quotient modules:", diagram3.generate_all_quotientmodules_bitmask())

    # Example 4: weird radical
    print("\n--- Example 4 ---")
    arrows4 = [
        Arrow(0, 1, "a"),
        Arrow(1, 2, "b"),
        Arrow(2, 3, "c"),
        Arrow(3, 4, "d"),
        Arrow(0, 5, "e"),
        Arrow(5, 4, "f")
    ]
    diagram4 = ModuleDiagram(arrows4)
    diagram4.draw_radical_layers()

    print("Nodes:", diagram4.nodes())
    print("Radical layers:", diagram4.node_to_radical_layers())
    print("All submodules:", diagram4.generate_all_submodules_bitmask())
    print("All quotient modules:", diagram4.generate_all_quotientmodules_bitmask())

    # Example 5: product of modules
    print("\n--- Example 5: product of modules ---")
    arrows5 = [
        Arrow(0, 1, "a1"),
        Arrow(1, 2, "b1"),
        Arrow(2, 3, "c1"),
        Arrow(4, 5, "a2"),
        Arrow(5, 6, "b2"),
        Arrow(6, 7, "c2")
    ]
    diagram5 = ModuleDiagram(arrows5)
    diagram5.draw_radical_layers()

    print("Nodes:", diagram5.nodes())
    print("Radical layers:", diagram5.node_to_radical_layers())
    print("All submodules:", diagram5.generate_all_submodules_bitmask())
    print("All quotient modules:", diagram5.generate_all_quotientmodules_bitmask())
    
