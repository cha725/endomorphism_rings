from quiver_algebra import MonomialQuiverAlgebra, Arrow, Path
import networkx as nx
from modulediagram import ModuleDiagram
import matplotlib.pyplot as plt

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
        if top_vertex not in algebra.vertices():
            raise ValueError(f"Invalid vertex. Vertex {top_vertex} must be a quiver vertex {algebra.vertices()}")
        self.algebra = algebra
        self.top_vertex = top_vertex
        self._construct_mod_diagram()
        super().__init__(arrows=self.arrows, 
                         vertex_simples=self.vertex_simples, 
                         isolated_vertices=self.isolated_vertices)

    def _construct_mod_diagram(self):
        self.vertex_simples = {}
        self.arrows = []
        self.isolated_vertices = None
        path_to_vertex_id = {}
        for idx, connection in enumerate(self.algebra.paths(with_connections=True)[self.top_vertex]):
            path, arrow, concatenation = connection
            path_to_vertex_id[concatenation] = idx
            self.vertex_simples[idx] = concatenation.target()
            if arrow is not None:
                source = path_to_vertex_id[path]
                target = path_to_vertex_id[concatenation]
                label = arrow
                self.arrows.append(Arrow(source,target,label))
        if not self.arrows and self.vertex_simples:
            self.isolated_vertices = list(range(len(self.vertex_simples))) 



if __name__ == "__main__":

    qa = MonomialQuiverAlgebra([Arrow(0,1), Arrow(1,2), Arrow(2,3)])
    for vertex in qa.vertices():
        print(ProjectiveModuleDiagram(qa, vertex))

    qa2 = MonomialQuiverAlgebra(
    [Arrow(0,1),Arrow(1,2),Arrow(2,0)],
    [Path((Arrow(0,1),Arrow(1,2)))]
    )
    for vertex in qa2.vertices():
        print(ProjectiveModuleDiagram(qa2, vertex))


