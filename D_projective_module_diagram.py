from A_bitmask_subgraph import Vertex, Arrow
from B_quiver_algebra import Path, MonomialQuiverAlgebra
from C_module_diagram import ModuleDiagram


class PathModuleDiagram(ModuleDiagram):
    """
    Represents a module diagram whose vertices correspond to paths in a path algebra.

    Attributes:
        - vertex_list (list[Vertex]):
            List of vertices, their labels are assumed to be paths in the path algebra.
        - arrow_list (list[Arrow]):
            List of arrows, their labels are assumed to be arrows in the path algebra.
    """
    def __init__(self,
                 vertex_list: list[Vertex],
                 arrow_list: list[Arrow] | None = None):
        self.vertex_list = vertex_list
        self.arrow_list = arrow_list or []
        self.reduce_labels()
        super().__init__(self.vertex_list,self.arrow_list)

    def strip_once(self) -> bool:
        """
        Each vertex label is a Path object. If each of these paths starts
        with the same arrow remove that arrow from each vertex label.

        Returns:
            - True if an arrow was stripped from the front.
            - False otherwise.
        """
        if not self.vertex_list:
            return False
        self.vertex_list.sort(key=lambda v:len(v.label))
        smallest_path: Path = self.vertex_list[0].label
        if len(smallest_path) == 0:
            # Then smallest path is a stationary path and cannot be reduced
            return False
        arrow_to_remove = smallest_path.first_arrow()
        old_vertex_to_new_vertex = {}
        for old_vertex in self.vertex_list:
            label: Path = old_vertex.label
            first_arrow = label.first_arrow()
            if arrow_to_remove != first_arrow:
                return False
            if len(label) == 1:
                new_label = Path(label.target())
            else:
                new_label = label.truncate(1)
            composition_factor = old_vertex.composition_factor
            radical_layer = old_vertex.radical_layer
            socle_layer = old_vertex.socle_layer
            new_vertex = Vertex(
                new_label,
                composition_factor,
                radical_layer,
                socle_layer
                )
            old_vertex_to_new_vertex[old_vertex] = new_vertex
        self.vertex_list = list(old_vertex_to_new_vertex.values())
        old_arrow_to_new_arrow = {}
        for arrow in self.arrow_list:
            new_source = old_vertex_to_new_vertex[arrow.source]
            new_target = old_vertex_to_new_vertex[arrow.target]
            old_arrow_to_new_arrow[arrow] = Arrow(new_source,new_target,arrow.label)
        self.arrow_list = list(old_arrow_to_new_arrow.values())
        return True

    def reduce_labels(self):
        """ 
        Returns a set of vertices whose labels have been fully reduced.
            i.e. They do not have a common first arrow.
        """
        while True:
            changed = self.strip_once()
            if not changed:
                break

    def __eq__(self, other: "PathModuleDiagram"):
        """ 
        The arrows of a PathModuleDiagram are determined by the vertex labels.
        Fix two vertices v with label p and u with label q. There is an arrow
        p -> q if and only if there q = pa for some arrow a.
        """
        return self.vertex_list == other.vertex_list and self.arrow_list == other.arrow_list


class ProjectiveDiagram(PathModuleDiagram):
    """
    A module diagram representing an indecomposable projective module
    of a MonomialQuiverAlgebra.

    Attributes:
        - algebra (MonomialQuiverAlgebra):
            The algebra over which to compute the projective module.
        - top_vertex (Vertex):
            The vertex of the quiver corresponding to the top of the projective.
    """
    def __init__(self,
                 algebra: MonomialQuiverAlgebra,
                 top_vertex: Vertex):
        if top_vertex not in algebra.vertices:
            raise ValueError(f"Invalid vertex. Vertex {top_vertex} must be a quiver vertex {algebra.vertices}")
        self.algebra = algebra
        self.top_vertex = top_vertex
        vertex_list, arrow_list = self._create_module_diagram()
        super().__init__(vertex_list, arrow_list)

    def _create_module_diagram(self) -> tuple[list[Vertex], list[Arrow]]:
        """
        Construct the vertex and arrow list for the projective module diagram.

        The list of vertices in the projective module diagram is the list of 
        paths in the quiver algebra starting at the top vertex. The vertex label
        is the path label.

        There is an arrow from a vertex p to a vertex q if there exists an arrow
        a such that pa = q. The label of this arrow is a.

        Returns:
            - tuple(list[Vertex], list[Arrow]):
                List of vertices labelled by paths in the algebra.
                List of arrows labelled by the arrows in the algebra.
        Note:
            - No two vertices will have the same label.
            - More than one arrow can have the same label.
        """
        paths = self.algebra.dfs_paths_from_vertex(self.top_vertex)
        assert paths is not None
        paths_to_vertex = {
            path : Vertex(label = path, composition_factor=path.target()) 
            for path in paths
        }
        vertex_list: list[Vertex] = list(paths_to_vertex.values())
        arrow_list: list[Arrow] = []
        for path in paths:
            if path.arrows:
                final_arrow = path.arrows[-1]
                previous_path = path.truncate(0,len(path)-1)
                source = paths_to_vertex[previous_path]
                target = paths_to_vertex[path]
                label = final_arrow.label
                arrow_list.append(Arrow(source, target, label))
        return (vertex_list, arrow_list)
    
    def nice_print(self):
        s = f"\nVertices = {[v.label for v in self.vertex_list]}"
        s += f"\nArrows = {[a.label for a in self.arrow_list]}"
        return s

    def __repr__(self):
        return f"Projective(algebra={self.algebra}, top vertex={self.top_vertex})"

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

    class Examples:
        """
        Class to store examples.
        """
        def __init__(self,
                        examples: dict[str,MonomialQuiverAlgebra]):
            self.examples = examples

        def add(self, example: tuple[str, MonomialQuiverAlgebra]):
            self.examples[example[0]] = example[1]

        def run(self, draw: bool = False):
            for name, quiver in self.examples.items():
                print(f"\n=== Example: {name} ===")
                print(f"\n Vertices of quiver = {[v.label for v in quiver.vertices]}")
                print(f" Arrows of quiver = {quiver.arrows}")
                print(f" Relations = {quiver.relations}")
                print(f"\n-- Projectives of quiver --")
                for v in quiver.vertices:
                    print(f"\n=== Projective at {v} ===")
                    proj = ProjectiveDiagram(quiver, v)
                    print(f"\n {proj.nice_print()}")
                    print(f"\n Radical layers of proj: {[[v.label for v in layer] for layer in proj.radical_layer_list]}")
                    print(f"\n Radical subgraphs:")
                    for r in proj.radical_submodules:
                        print([v.label for v in r.vertex_list])                
                    if draw:
                        proj.draw_radical_layers

    examples = Examples({})

    examples.add(("Type A no relations",
                  MonomialQuiverAlgebra(arrows=[Arrow(Vertex(0),Vertex(1),"a"),
                                                Arrow(Vertex(1),Vertex(2),"b"),
                                                Arrow(Vertex(2),Vertex(3),"c")],
                                        relations = [])))
    
    examples.add(("Type A rad2 relations",
                  MonomialQuiverAlgebra(arrows=[Arrow(Vertex(0),Vertex(1),"a"),
                                                Arrow(Vertex(1),Vertex(2),"b"),
                                                Arrow(Vertex(2),Vertex(3),"c")],
                                        relations = [Path((Arrow(Vertex(0),Vertex(1),"a"),
                                                           Arrow(Vertex(1),Vertex(2),"b"))),
                                                     Path((Arrow(Vertex(1),Vertex(2),"b"),
                                                           Arrow(Vertex(2),Vertex(3),"c")))])))

    examples.add(("Three cyclic",
                 MonomialQuiverAlgebra(arrows=[Arrow(Vertex(0),Vertex(1),"a"),
                                               Arrow(Vertex(1),Vertex(2),"b"),
                                               Arrow(Vertex(2),Vertex(0),"c")],
                                       relations=[Path((Arrow(Vertex(0),Vertex(1),"a"),
                                                        Arrow(Vertex(1),Vertex(2),"b")))])))

    examples.add(("Bigger",
                 MonomialQuiverAlgebra( arrows=[Arrow(Vertex(1),Vertex(2),"a"),
                                                Arrow(Vertex(2),Vertex(3),"b"),
                                                Arrow(Vertex(1),Vertex(3),"c"),
                                                Arrow(Vertex(3),Vertex(1),"d"),
                                                Arrow(Vertex(2),Vertex(2),"e")],
                                        relations=[Path((Arrow(Vertex(1),Vertex(2),"a"),
                                                         Arrow(Vertex(2),Vertex(2),"e")))],
                                        max_radical_length=3)))

    examples.run()
