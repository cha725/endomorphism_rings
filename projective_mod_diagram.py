from quiver_algebra import MonomialQuiverAlgebra
from modulediagram import ModuleDiagram

class ProjectiveModuleDiagram(ModuleDiagram):
    def __init__(self,
                 algebra : MonomialQuiverAlgebra,
                 top_vertex : int):
        if top_vertex not in algebra.vertices():
            raise ValueError(f"Invalid vertex. Vertex {top_vertex} must be a quiver vertex {algebra.vertices()}")
        self.algebra = algebra
        self.top_vertex = top_vertex
        self.arrows = self.algebra.paths_from_fixed_pt(self.top_vertex,self.algebra.max_radical_length)
        # TODO: need to track the composition from the quiver algebra vertices.
        # TODO: think of a a better word than composition.
        super().__init__(self.arrows)

