from quiver_algebra import MonomialQuiverAlgebra, Arrow, Path
from modulediagram import ModuleDiagram

class ProjectiveModuleDiagram(ModuleDiagram):
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
        vertex_to_idx = {}
        composition_factors = []
        vertex_labels = []
        arrows = []

        dfs_connections = self.algebra.paths(with_connections=True)[self.top_vertex]

        for path, arrow, new_path in dfs_connections:
            print(f"path {path}, new path {new_path}")

            if new_path not in vertex_to_idx:
                idx = len(vertex_to_idx) # take next index available
                vertex_to_idx[new_path] = idx # assign new path this index
                composition_factors.append(new_path.target()) # assign simple factor to this index
                vertex_labels.append(new_path) # 

            if arrow is not None:
                source = vertex_to_idx[path]
                target = vertex_to_idx[new_path]
                arrows.append(Arrow(source,target,arrow.label if hasattr(arrow, 'label') else arrow))

        self.composition_factors = tuple(composition_factors)
        self.vertex_labels = tuple(vertex_labels)
        self.arrows = tuple(arrows)





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

    def run(self):
        for name, quiver in self.examples.items():
            print(f"\n=== Example: {name} ===")
            print(f"\n Vertices = {quiver.vertices}")
            print(f" Arrows = {quiver.arrows}")
            print(f" Relations = {quiver.relations}")
            print(f"\n-- Projectives of quiver --")
            for vertex in quiver.vertices:
                projective = ProjectiveModuleDiagram(quiver, vertex)
                print(f"\n- Projective {projective} -")
                projective.draw_radical_layers
            for submod in projective.generate_all_submodules:
                print(f"submod {submod}")

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
