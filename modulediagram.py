import matplotlib.pyplot as plt
import networkx as nx
import random
import time

from functools import cached_property
from types import MappingProxyType

from bitmask_subgraph import Vertex, Arrow, BitmaskSubgraph

class ModuleDiagram:
    """
    Represents a directed graph modelling a module.

    - Each vertex corresponds to a composition factor of the module.
    - Each arrow represents the action of a ring element on the composition factor.

    Provides functionality to compute radical layers, submodules, and quotients.

    Note: the directed graph is acyclic, but the undirected graph may not be acyclic.

    Parameters:
        vertex_list (list[Vertex]): List of vertices of the module diagram.
        arrow_list (list(Arrow,...)): List of arrows between the vertices of the module diagram.
    """
    def __init__(self, 
                 vertex_list: list[Vertex] | None = None,
                 arrow_list: list[Arrow] | None = None):
        
        # Compute vertices
        if vertex_list:
            self.vertex_list = vertex_list
        elif arrow_list:
            unordered = set(a.source for a in arrow_list) | set(a.target for a in arrow_list)
            try:
                self.vertex_list = sorted(unordered, key=lambda v : v.label)
            except(TypeError, AttributeError):
                self.vertex_list = list(unordered)
        else:
            raise ValueError("Either vertex_list or arrow_list must be provided.")
        
        self.vertex_labels = [v.label for v in self.vertex_list]
        self.num_vertices = len(self.vertex_list)

        self.vertex_to_index: MappingProxyType[Vertex, int] = MappingProxyType({v: idx for idx, v in enumerate(self.vertex_list)})
        self.index_to_vertex: MappingProxyType[int, Vertex] = MappingProxyType({idx: v for idx, v in enumerate(self.vertex_list)})

        # Compute arrows
        
        if arrow_list:
            for a in arrow_list:
                if a.source not in self.vertex_list:
                    raise ValueError(f"Arrow sources must be in vertex list. Got arrow {a} for vertex list {self.vertex_list}.")
                if a.target not in self.vertex_list:
                    raise ValueError(f"Arrow targets must be in vertex list. Got arrow {a} for vertex list {self.vertex_list}.")
        self.arrow_list = arrow_list or []
        self.arrow_labels = [a.label for a in self.arrow_list]
        self.num_arrows = len(self.arrow_list)

        self.bitmask = BitmaskSubgraph(self.vertex_list, self.arrow_list)
    
    def sources(self) -> list[Vertex]:
        """ Return list of vertices with no incoming arrows. """
        return self.bitmask.sources()

    def sinks(self) -> list[Vertex]:
        """ Return list of vertices with no outgoing arrows. """
        return self.bitmask.sinks()

    @cached_property
    def radical_layer_list(self) -> list[list[Vertex]]:
        """ Return list of vertices grouped into radical layers. """
        return self.bitmask.radical_layers()

    @cached_property
    def radical_layer_of(self) -> dict[Vertex, int]:
        """ Return dictionary assigning to each vertex its radical layer. """
        radical_labels = {}
        for layer_idx, layer in enumerate(self.radical_layer_list):
            for v in layer:
                radical_labels[v] = layer_idx
        return radical_labels

    
    def loewy_length(self) -> int:
        """
        Return Loewy length of module diagram.
        Loewy length is the smallest n such that the nth radical layer is trivial.
        """
        return len(self.radical_layer_list)
    
    @cached_property
    def socle_layers(self) -> list[list[Vertex]]:
        """ Return list of vertices grouped into socle layers. """
        return self.bitmask.socle_layers()

    @cached_property
    def socle_layer_of(self) -> dict[Vertex, int]:
        """ Return dictionary assigning to each vertex its socle layer. """
        socle_labels = {}
        for layer_idx, layer in enumerate(self.socle_layers):
            for v in layer:
                socle_labels[v] = layer_idx
        return socle_labels
    @cached_property
    def radical_layer_to_vertices(self) -> list[list[int]]:
        """
        Returns:
            list[list[int]]: A list entry i = all vertices in radical layer i.
        """
        r_labels = self.radical_labels
        r_layers = [[] for _ in range(self.num_radical_layers())]
        for v, lay in enumerate(r_labels):
            r_layers[lay].append(v)
        return r_layers

    def _create_radical_layer_graph(self):
        """
        Create graph with vertex labels: composition factor, vertex label, radical layer
        """
        G = nx.MultiDiGraph()
        for idx, v in enumerate(self.composition_factors):
            if self.vertex_labels:
                v_label = self.vertex_labels[idx]
            r_label = self.radical_labels[idx]
            G.add_node(idx, comp_factor=v, label=v_label, rad_layer=r_label)
        for arrow in self.arrows:
            G.add_edge(arrow.source, arrow.target, label=arrow.label)
        return G

    @cached_property
    def draw_radical_layers(self):
        """
        Draw the nxgraph lined up with radical layers.
        """
        pos = {}
        for layer, vertices in enumerate(self.radical_layer_to_vertices):
            n = len(vertices)
            if n == 1:
                x_coord = [0.5]
            else:
                x_coord = [i / (n - 1) for i in range(n)]
            for x, node in zip(x_coord, vertices):
                pos[node] = (x, -layer)

        G = self._create_radical_layer_graph()
        nx.draw(G,
                pos=pos,
                node_size=1000, 
                node_color='lightblue',
                font_weight='bold')
        nx.draw_networkx_labels(G, 
                                pos, 
                                labels={n: G.nodes[n]['comp_factor'] for n in G.nodes},
                                font_size=10, 
                                font_weight='bold')  
        edge_labels = {
            (u, v, k): d.get("label")
            for u, v, k, d in self.graph.edges(keys=True, data=True)
            if d.get("label") is not None
        }

        nx.draw_networkx_edge_labels(
            self.graph, pos,
            edge_labels=edge_labels,
            font_size=9,
            label_pos=0.5,
            rotate=False
        )      
        plt.show()



    def compute_isomorphism(self, other: "ModuleDiagram") -> dict | None:
        G = self._create_radical_layer_graph
        H = other._create_radical_layer_graph
        matcher = nx.algorithms.isomorphism.DiGraphMatcher(
            G,
            H,
            node_match = lambda a, b: a.get("comp_factor") == b.get("comp_factor"),
            edge_match = lambda a, b: a.get("label") == b.get("label")
        )
        if not matcher.is_isomorphic():
            return None
        return dict(matcher.mapping)

    def __eq__(self, other : "ModuleDiagram") -> bool:
        return self.compute_isomorphism(other) is not None
    
    def __repr__(self):
        return f"MD(vertices = {self.vertex_labels}, arrows = {self.arrow_list})"

class ModuleSubDiagram(ModuleDiagram):
    def __init__(self,
                 parent : ModuleDiagram,
                 vertex_sublist : list[Vertex]):
        for v in vertex_sublist:
            if v not in parent.vertex_list:
                raise ValueError(f"{v} is not a vertex of parent module, {parent.vertex_list}")
        self.parent = parent

        if parent.arrow_list:
            arrow_sublist = [
                arr for arr in parent.arrow_list
                    if arr.source in vertex_sublist and arr.target in vertex_sublist
            ]
        else:
            arrow_sublist = []

        super().__init__(
            vertex_sublist,
            arrow_sublist
        )

    def compute_isomorphism(self, other: "ModuleSubDiagram") -> dict[Vertex, Vertex] | None:
        """
        Returns a mapping using parent vertex indices if isomorphism is found.
        Or None if an isomorphism is not found.
        """
        G = self._create_radical_layer_graph
        H = other._create_radical_layer_graph

        matcher = nx.algorithms.isomorphism.DiGraphMatcher(
            G,
            H,
            node_match = lambda a, b: a.get("comp_factor") == b.get("comp_factor"),
            edge_match = lambda a, b: a.get("label") == b.get("label")
        )

        if not matcher.is_isomorphic():
            return None

        return matcher.mapping

    def is_submodule(self) -> bool:
        dim = self.num_vertices
        return self.vertex_list in self.parent.all_submodules[dim]
    
    def is_quotient(self) -> bool:
        dim = self.num_vertices
        return self.vertex_list in self.parent.all_quotients[dim]
    

    def add(self, 
            name: str, 
            composition_factors: tuple[int,...],
            vertex_labels: Optional[tuple[str,...]] = None,
            arrows: Optional[tuple[Arrow,...]] = None):
        self.examples[name] = { "composition_factors" : composition_factors,
                                "vertex_labels" : vertex_labels,
                                "arrows" : arrows}

    def add_random_diagram(self, name: str, num_vertices : int, edge_prob: float):
        """
        Add a random directed tree with num_vertices vertices and edge probability edge_prob.
        """
        composition_factors = tuple([0]*num_vertices)
        vertex_labels = tuple(['']*num_vertices)
        arrows = []
        for i in range(num_vertices):
            for j in range(i+1,num_vertices):
                if random.random() < edge_prob:
                    arrows.append(Arrow(i, j))
        arrows = tuple(arrows)
        self.add(name, composition_factors, vertex_labels, arrows)

    def run(self, verbose : bool = False, draw: bool = False):
        times = []
        for name, data in self.examples.items():
            print(f"\n=== Example: {name} ===")

            start_time = time.time()
            diagram = ModuleDiagram(composition_factors=data["composition_factors"],
                                    vertex_labels=data["vertex_labels"],
                                    arrows=data["arrows"])
            init_time = time.time() - start_time

            start_time = time.time()
            submods = diagram.generate_all_submodules
            quotients = diagram.generate_all_quotients
            subgraph_time = time.time() - start_time

            if draw:
                diagram.draw_radical_layers

            if verbose:
                print("\nVertices:", diagram.vertices)
                print("\nArrows:", diagram.arrows)
                print("\nSources:", diagram.sources)
                print("\nRadical layers:", diagram.radical_layer_to_vertices)
                print("\nAll submodules:", submods)
                print("\nAll quotient modules:", quotients)
        
            times.append((init_time, subgraph_time))

        return times

# TODO: make the print out prettier.

if __name__ == "__main__":

    examples = Examples({})

    for n in range(1):
        examples.add_random_diagram(f"Example {n}", 10, 0.4)

    times = examples.run(draw=True, verbose=True)

    for idx, t in enumerate(times):
        print(f"\n=== Example {idx} ===")
        print(f"Initialisation time: {t[0]:.4f}")
        print(f"Subgraph time: {t[1]:.4f}")




    # examples.add("Example 1: A4",
    #                     arrows=[Arrow(0,1,"a"), Arrow(1,2,"b"), Arrow(2,3,"c")],
    #                     vertex_simples={0:0, 1:1, 2:0, 3:1})

    # examples.add("Example 2: Diamond",
    #                     arrows=[Arrow(0,1,"a"), Arrow(0,2,"b"), Arrow(1,3,"c"), Arrow(2,3,"d")])

    # examples.add("Example 3: Y structure",
    #                     arrows=[Arrow(0,2,"a"), Arrow(1,2,"b"), Arrow(2,3,"c"), Arrow(3,4,"d")])

    # examples.add("Example 4: weird radical",
    #                     arrows=[Arrow(0,1,"a"), Arrow(1,2,"b"), Arrow(2,3,"c"), Arrow(3,4,"d"),
    #                             Arrow(0,5,"e"), Arrow(5,4,"f")])

    # examples.add("Example 5: product of modules",
    #                     arrows=[Arrow(0,1,"a1"), Arrow(1,2,"b1"), Arrow(2,3,"c1"),
    #                             Arrow(4,5,"a2"), Arrow(5,6,"b2"), Arrow(6,7,"c2")])

    # examples.add("Example 6: Algebra and dual",
    #                     arrows=[Arrow(0,1,"a1"), Arrow(0,2,"b1"), Arrow(4,5,"a2"), Arrow(6,7,"b2")],
    #                     isolated_vertices=[3])
    

