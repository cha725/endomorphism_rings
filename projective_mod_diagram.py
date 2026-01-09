from bitmask_subgraph import Vertex, Arrow
from quiver_algebra import MonomialQuiverAlgebra
from modulediagram import ModuleDiagram




class ProjectiveDiagram(ModuleDiagram):
    """
    A module diagram representating an indecomposable projective module
    of a MonomialQuiverAlgebra.

    Attributes:
        algebra (MonomialQuiverAlgebra):
            The algebra over which to compute the projective module.
        top_vertex (int):
            The vertex of the quiver corresponding to the top of the projective.
    """
    def __init__(self,
                 algebra : MonomialQuiverAlgebra,
                 top_vertex : int):
        if top_vertex not in algebra.vertices:
            raise ValueError(f"Invalid vertex. Vertex {top_vertex} must be a quiver vertex {algebra.vertices}")
        self.algebra = algebra
        self.top_vertex = top_vertex
        vertex_list, arrow_list = self._create_module_diagram()
        super().__init__(vertex_list, arrow_list)

    def _create_module_diagram(self) -> tuple[list[Vertex], list[Arrow]]:
        paths, connections = self.algebra.dfs_paths_from_vertex(self.top_vertex)
        paths_to_vertex = {path : Vertex(label = path.fancy_label(), composition_factor=path.target()) for path in paths}
        vertex_list = list(paths_to_vertex.values())
        arrow_list = []
        for connection in connections:
            source = paths_to_vertex[connection.source]
            target = paths_to_vertex[connection.target]
            label = connection.label
            arrow_list.append(Arrow(source, target, label))
        return (vertex_list, arrow_list)
    

    def nice_print(self):
        s = f"\nVertices = {self.vertex_list}"
        s += f"\nArrows = {[a.label for a in self.arrow_list]}"
        return s

    def __repr__(self):
        return f"Projective(algebra={self.algebra}, top vertex={self.top_vertex})"


# TODO: Is this standard practice?

class Examples:
    """
    Class to store examples.
    """
    def __init__(self,
                    examples : dict[str,MonomialQuiverAlgebra]):
        self.examples = examples

    def add(self, example : tuple[str, MonomialQuiverAlgebra]):
        self.examples[example[0]] = example[1]

    def run(self, draw: bool = False):
        for name, quiver in self.examples.items():
            print(f"\n=== Example: {name} ===")
            print(f"\n Vertices = {quiver.vertices}")
            print(f" Arrows = {quiver.arrows}")
            print(f" Relations = {quiver.relations}")
            print(f"\n-- Projectives of quiver --")
            for v in quiver.vertices:
                proj = ProjectiveDiagram(quiver, v)
                print(proj)
                if draw:
                    proj.draw_radical_layers
                print(f"Paths = {proj.paths}")
                print(f"Submodules:")
                for submod in proj.find_all_submodules():
                    print(f"{submod}")

if __name__ == "__main__":

    examples = Examples({})

    examples.add(("Type A no relations",
                  MonomialQuiverAlgebra(arrows=[Arrow(0,1,"a"),Arrow(1,2,"b"),Arrow(2,3,"c")],
                                        relations = [])))
    
    examples.add(("Type A rad2 relations",
                  MonomialQuiverAlgebra(arrows=[Arrow(0,1,"a"),Arrow(1,2,"b"),Arrow(2,3,"c")],
                                        relations = [Path((Arrow(0,1,"a"),Arrow(1,2,"b"))),
                                                     Path((Arrow(1,2,"b"),Arrow(2,3,"c")))])))

    examples.add(("Three cyclic",
                 MonomialQuiverAlgebra(arrows=[Arrow(0,1,"a"),Arrow(1,2,"b"),Arrow(2,0,"c")],
                                       relations=[Path((Arrow(0,1,"a"),Arrow(1,2,"b")))])))

    examples.add(("Bigger",
                 MonomialQuiverAlgebra( arrows=[Arrow(1,2,"a"), Arrow(2,3,"b"), Arrow(1,3,"c"),
                                                Arrow(3,1,"d"), Arrow(2,2,"e")],
                                        relations=[Path((Arrow(1,2,"a"), Arrow(2,2,"e")))],
                                        max_radical_length=3)))

    examples.run()
