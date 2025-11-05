from quiver_algebra import MonomialQuiverAlgebra, Arrow, Path
import networkx as nx
from modulediagram import ModuleDiagram
import matplotlib.pyplot as plt

class ProjectiveModuleDiagram(ModuleDiagram):
    def __init__(self,
                 algebra : MonomialQuiverAlgebra,
                 top_vertex : int):
        if top_vertex not in algebra.vertices():
            raise ValueError(f"Invalid vertex. Vertex {top_vertex} must be a quiver vertex {algebra.vertices()}")
        self.algebra = algebra
        self.top_vertex = top_vertex
        self.vertex_simples = {}
        self.arrows = []
        path_to_vertex_id = {}
        for idx, connection in enumerate(self.algebra.paths(with_connections=True)[top_vertex]):
            if len(connection) < 3:
                path_to_vertex_id[connection] = idx
                self.vertex_simples[idx] = connection.target()
            else:
                path_to_vertex_id[connection[2]] = idx
                self.vertex_simples[idx] = connection[2].target()
                self.arrows.append(Arrow(path_to_vertex_id[connection[0]],target=path_to_vertex_id[connection[2]],label=connection[1]))
        super().__init__(arrows=self.arrows, vertex_simples=self.vertex_simples)
        # TODO: Not picking up the projective modules that are 1 dimensional
        # Thinks there are no simples in them.


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


