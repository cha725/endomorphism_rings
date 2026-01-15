from A_bitmask_subgraph import Vertex, Arrow
from B_quiver_algebra import Path, MonomialQuiverAlgebra
from C_module_diagram import ModuleDiagram




class ProjectiveDiagram(ModuleDiagram):
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
            path : Vertex(label = path.label(), composition_factor=path.target()) 
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
