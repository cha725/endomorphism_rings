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

    # @cached_property
    # def radical_layers(self) -> list[list[Vertex]]:
    #     """ Return list of vertices grouped into radical layers. """
    #     rad_layers_idx = self._radical_layers_idx
    #     return [[self.index_to_vertex[idx] for idx in layer] for layer in rad_layers_idx]

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
    
    ## Submodules and quotient modules ##

    @cached_property
    def all_submodules(self) -> list[list[list[Vertex]]]:
        """ Return all submodules (successor-closed subsets) indexed by size. """
        return self.bitmask.compute_succ_closed_subsets()
    
    @cached_property
    def all_quotients(self) -> list[list[list[Vertex]]]:
        """ Return all quotient modules (predeccesor-closed subsets) indexed by size. """
        return self.bitmask.compute_pred_closed_subsets()

    
    # def is_submodule(self, vertex_sublist: list[Vertex]) -> bool:
    #     """
    #     Check if given vertex sublist forms a submodule of self.
    #     A submodule corresponds to a subset of vertices that is successor-closed:
    #         i.e. the successors of every vertex in the subset are also in the subset.
    #     """
    #     dim = len(vertex_sublist)
    #     return vertex_sublist in self.all_submodules[dim]
    
    # def is_quotient(self, vertex_sublist: list[Vertex]) -> bool:
    #     """
    #     Check if given vertex sublist forms a quotient module of self.
    #     A quotient module corresponds to a subset of vertices that is predecessor-closed:
    #         i.e. the predecessors of every vertex in the subset are also in the subset.
    #     """
    #     dim = len(vertex_sublist)
    #     return vertex_sublist in self.all_quotients[dim]
    
    @cached_property
    def radical_submodules(self) -> list["ModuleSubDiagram"]:
        """ Returns list of radical power submodules. """
        vertex_lists = self.bitmask.compute_radical_subgraphs()
        return [ModuleSubDiagram(self, v_list) for v_list in vertex_lists]
    
    def indecomposable_summands(self) -> "list[ModuleSubDiagram]":
        """ Returns list of indecomposable summands. """
        subsets = self.bitmask.connected_components
        return [ModuleSubDiagram(self, v) for v in subsets]



    ## Drawing the module diagram ##
    @cached_property
    def _create_radical_layer_graph(self):
        """ 
        Construct a networkx directed graph of module.
        Node labels: 
            - label = vertex label of vertex (if given).
            - comp_factor = composition factor represented by vertex.
            - rad_layer = radical layer of vertex.
        """
        G = nx.MultiDiGraph()
        for v in self.vertex_list:
            G.add_node(v, comp_factor = v.composition_factor, label = v.label, rad_layer = v.radical_layer)
        for arrow in self.arrow_list:
            G.add_edge(arrow.source, arrow.target, label=arrow.label)
        return G

    @cached_property
    def draw_radical_layers(self):
        """ Draw the module diagram with nodes aligned by radical layer. """
        pos = {}
        for layer_idx, vertices in enumerate(self.radical_layer_list):
            n = len(vertices)
            if n == 1:
                x_coord = [0.5]
            else:
                x_coord = [i / (n - 1) for i in range(n)]
            for x, v in zip(x_coord, vertices):
                pos[v] = (x,-layer_idx)

        G = self._create_radical_layer_graph
        nx.draw(G,
                pos=pos,
                node_size=1000, 
                node_color='lightblue',
                font_weight='bold')
        nx.draw_networkx_labels(G, 
                                pos, 
                                labels={n: G.nodes[n]["comp_factor"] for n in G.nodes},
                                font_size=10, 
                                font_weight='bold')  
        edge_labels = {
            (u, v, k): d.get("label")
            for u, v, k, d in G.edges(keys=True, data=True)
            if d.get("label") is not None
        }

        nx.draw_networkx_edge_labels(
            G, pos,
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
    
    __hash__ = object.__hash__


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
    
    def indecomposable_summands(self) -> list[ModuleSubDiagram]:
        """ Returns list of indecomposable summands as subdiagrams of the parent module. """
        ind_summands = super().indecomposable_summands()
        change_parent = []
        for summand in ind_summands:
            change_parent.append(ModuleSubDiagram(self.parent, summand.vertex_list))
        return change_parent

    def is_submodule(self) -> bool:
        dim = self.num_vertices
        return self.vertex_list in self.parent.all_submodules[dim]
    
    def is_quotient(self) -> bool:
        dim = self.num_vertices
        return self.vertex_list in self.parent.all_quotients[dim]
    
    def __repr__(self):
        return f"MSubD(parent={self.parent}, vertices={self.vertex_labels})"
    
    __hash__ = object.__hash__



###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################




if __name__ == "__main__":


    class MDExamples:
        """
        Class to store examples of module diagrams.
        """
        def __init__(self, 
                    examples: dict[str, dict]):
            self.examples = examples

        def add(self, 
                name: str,
                vertex_list: list[Vertex] | None = None,
                arrow_list: list[Arrow] | None = None):
            self.examples[name] = { 
                "vertex_list" : vertex_list,
                "arrow_list" : arrow_list
            }

        def add_random_diagram(self, name: str, num_vertices : int, edge_prob: float):
            """
            Add a random directed tree with num_vertices vertices and edge probability edge_prob.
            """
            vertex_list = [Vertex(chr(i+97)) for i in range(num_vertices)]
            arrow_list = []
            for i in range(num_vertices):
                for j in range(i+1,num_vertices):
                    if random.random() < edge_prob:
                        arrow_list.append(Arrow(Vertex(chr(i+97)), Vertex(chr(j+97))))
            self.add(name, vertex_list, arrow_list)

        def run(self, verbose : bool = False, draw: bool = False):
            times = []
            for name, data in self.examples.items():
                print(f"\n=== Example: {name} ===")

                start_time = time.time()
                diagram = ModuleDiagram(
                    vertex_list=data["vertex_list"],
                    arrow_list=data["arrow_list"]
                )
                init_time = time.time() - start_time
                print(f"Time to initialise object: {init_time:.4f} seconds")

                start_time = time.time()
                decomposition = diagram.indecomposable_summands()
                decomp_time = time.time() - start_time
                print(f"Time to decompose object {decomp_time:.4f} seconds")

                start_time = time.time()
                rad_layers = diagram.radical_layer_list
                rad_submods = diagram.radical_submodules
                rad_time = time.time() - start_time
                print(f"Time to create radical subgraphs: {rad_time:.4f} seconds")

                start_time = time.time()
                subs = diagram.all_submodules
                quotients = diagram.all_quotients
                subgraph_time = time.time() - start_time
                print(f"Time to generate sub and quotient modules: {subgraph_time:.4f} seconds")

                if verbose:
                    print("\nVertices:", diagram.vertex_labels)
                    print("\nArrows:", diagram.arrow_list)
                    print("\nDecomposition:", [d.vertex_labels for d in decomposition])
                    print("\nRadical layers:", [[v.label for v in layer] for layer in rad_layers])
                    print("\nRadical submodules:", rad_submods)
                    print("\nSocle layers:", [[v.label for v in layer] for layer in diagram.socle_layers]) 
                    print("\nAll submodules:")
                    for s_list in subs:
                        for v_list in s_list:
                            print([v.label for v in v_list])
                    print("\nAll quotient modules:")
                    for q_list in quotients:
                        for v_list in q_list:
                            print([v.label for v in v_list])
            
                if draw:
                    diagram.draw_radical_layers
                times.append((init_time, subgraph_time))

            return times

    # TODO: make the print out prettier.

    class MSDExamples:
        """
        Class to store and run examples of ModuleSubDiagram instances
        based on a parent ModuleDiagram and vertex subsets.
        """
        def __init__(self, parent: ModuleDiagram):
            self.parent = parent
            self.examples: dict[str, list[Vertex]] = {}

        def add(self, name: str, vertices: list[Vertex]):
            """
            Add a subdiagram example using a subset of parent vertices.
            """
            if not set(vertices) <= set(self.parent.vertex_list):
                raise ValueError(f"Vertices {vertices} must be a subset of parent vertices {self.parent.vertex_list}")
            self.examples[name] = vertices

        def add_random_subdiagram(self, name: str):
            """
            Add a random directed tree with num_vertices vertices and edge probability edge_prob.
            """
            parent_verts = self.parent.vertex_list
            parent_num_verts = self.parent.num_vertices
            num_verts = random.choice(range(1, parent_num_verts))
            self.add(name, random.sample(parent_verts, num_verts))

        def run(self, verbose: bool = False):
            results = {}
            for name, vertices in self.examples.items():
                subdiagram = ModuleSubDiagram(self.parent, vertices)
                print("===========parent", self.parent)
                results[name] = subdiagram
                summands = subdiagram.indecomposable_summands()

                if verbose:
                    print(f"\n=== Subdiagram Example: {name} ===")
                    print(f"Parent diagram: {subdiagram.parent}")
                    print(f"Subdiagram vertices: {subdiagram.vertex_labels}")
                    print(f"Subdiagram arrows: {subdiagram.arrow_list}")
                    print(f"Idecomposable summands: {[s.vertex_labels for s in summands]}")
                    print(f"Which summands are submodules/quotients of parent module?")
                    for summand in summands:
                        print(f"Is {summand.vertex_list} a submodule of parent? {summand.is_submodule()}")
                        print(f"Is {summand.vertex_list} a quotient of parent? {summand.is_quotient()}")                   
                    
            return results


    # ### RUN EXAMPLES ###

    # ### ModuleDiagram Examples ###

    MD_examples = MDExamples({})

    for n in range(1):
       MD_examples.add_random_diagram(f"Example {n}", 7, 0.4)
    
    MD_examples.run(verbose=True, draw=True)

    # ### ModuleSubDiagram Examples ###

    # MSD_parent = ModuleDiagram(arrow_list=[Arrow(Vertex(0),Vertex(1)), Arrow(Vertex(1),Vertex(2))])

    MSD_parent_data = MD_examples.examples["Example 0"]
    MSD_parent_vertices = MSD_parent_data["vertex_list"]
    MSD_parent_arrows = MSD_parent_data["arrow_list"]
    MSD_parent = ModuleDiagram(MSD_parent_vertices, MSD_parent_arrows)

    MSD_examples = MSDExamples(MSD_parent)
    for n in range(10):
        MSD_examples.add_random_subdiagram(f"Example {n}")

    # MSD_examples.add("first two vertices", [Vertex(0),Vertex(1)])
    # MSD_examples.add("last two vertices", [Vertex(1),Vertex(2)])

    #MSD_examples.run(verbose=True)



    # # examples.add("Example 1: A4",
    # #                     arrows=(Arrow(0,1,"a"), Arrow(1,2,"b"), Arrow(2,3,"c")))

    # # examples.add("Example 2: Diamond",
    # #                     arrows=(Arrow(0,1,"a"), Arrow(0,2,"b"), Arrow(1,3,"c"), Arrow(2,3,"d")))

    # # examples.add("Example 3: Y structure",
    # #                     arrows=(Arrow(0,2,"a"), Arrow(1,2,"b"), Arrow(2,3,"c"), Arrow(3,4,"d")))

    # # examples.add("Example 4: weird radical",
    # #                     arrows=(Arrow(0,1,"a"), Arrow(1,2,"b"), Arrow(2,3,"c"), Arrow(3,4,"d"),
    # #                             Arrow(0,5,"e"), Arrow(5,4,"f")))

    # # examples.add("Example 5: product of modules",
    # #                     arrows=(Arrow(0,1,"a1"), Arrow(1,2,"b1"), Arrow(2,3,"c1"),
    # #                             Arrow(4,5,"a2"), Arrow(5,6,"b2"), Arrow(6,7,"c2")))

    # # examples.add("Example 6: Algebra and dual",
    # #                     arrows=(Arrow(0,1,"a1"), Arrow(0,2,"b1"), Arrow(4,5,"a2"), Arrow(6,7,"b2")),
    # #                     composition_factors=(0,1,2,3,4,5,6,7))
    
