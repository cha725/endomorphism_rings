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
        super().__init__(arrows=self.arrows, 
                         vertex_simples=self.vertex_simples, 
                         isolated_vertices=self.isolated_vertices,
                         vertex_labels=self.vertex_labels)

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
        self.vertex_simples = {} # node index to the simple factor
        self.vertex_labels = {} # node index to the path from dfs
        self.arrows = [] # possible arrows between nodes - really just arrows from the quiver
        self.isolated_vertices = []

        path_to_vertex_id = {}

        dfs_connections = self.algebra.paths(with_connections=True)[self.top_vertex]
        for path, arrow, new_path in dfs_connections:
            print(f"path {path}, new path {new_path}")

            if new_path not in path_to_vertex_id:
                idx = len(path_to_vertex_id) # take next index available
                path_to_vertex_id[new_path] = idx # assign new path this index
                self.vertex_simples[idx] = new_path.target() # assign simple factor to this index
                self.vertex_labels[idx] = new_path # 

            if arrow is not None:
                source = path_to_vertex_id[path]
                target = path_to_vertex_id[new_path]
                self.arrows.append(Arrow(source,target,arrow))

        if not self.arrows and self.vertex_simples:
            self.isolated_vertices = list(range(len(self.vertex_simples))) 



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
                #projective.draw_radical_layers
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