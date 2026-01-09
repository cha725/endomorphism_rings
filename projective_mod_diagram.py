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
        self._construct_mod_diagram()
        super().__init__(self.composition_factors,
                         self.vertex_labels,
                         self.arrows)

    def _construct_mod_diagram(self):
        """
        Build the projective module diagram for the given top vertex.
        
        Algebraically:
            - Each vertex corresponds to a simple composition factor of the module.
            - Each arrow represents the algebra action extending a path.
        
        Computationally:
            - Each vertex corresponds to a path in the DFS search of the quiver algebra.
            - Each arrow corresponds to a quiver arrow such that source_path + arrow = target_path.
        """
        path_to_vertex = {}
        composition_factors = []
        vertex_labels = []
        arrows = []
        paths = []

        dfs_connections = self.algebra.paths(with_connections=True)[self.top_vertex]

        for path, arrow, new_path in dfs_connections:

            if new_path not in path_to_vertex:
                idx = len(path_to_vertex) # take next index available
                path_to_vertex[new_path] = idx # assign new path this index
                composition_factors.append(new_path.target()) # assign simple factor to this index
                vertex_labels.append(new_path) # 
                paths.append(path)

            if arrow is not None:
                source = path_to_vertex[path]
                target = path_to_vertex[new_path]
                arrows.append(Arrow(source,target,arrow.label if hasattr(arrow, 'label') else arrow))

        self.path_to_vertex: dict[Path,int] = path_to_vertex
        self.composition_factors: tuple[int] = tuple(composition_factors)
        self.vertex_labels: tuple[str] = tuple(vertex_labels)
        self.arrows: tuple[Arrow] = tuple(arrows)
        self.paths: list[Path] = list(set(paths))

    def vertex_to_path(self) -> list[Path]:
        """
        Return list whose ith entry is the quiver path that vertex i is defined by.
        """
        vtp = [Path() for _ in range(len(self.path_to_vertex))]
        for v, path in enumerate(self.path_to_vertex):
            vtp[v] = path
        return vtp
    
    def find_all_submodules(self):
        submods = self.generate_all_submodules
        if submods is None:
            return None
        new_submods = []
        for submod in submods:
            new_submod = []
            for v in submod:
                new_submod.append(self.vertex_to_path()[v])
            new_submods.append(new_submod)
        return new_submods
    
    def find_all_submods_up_to_iso(self):
        pass

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
