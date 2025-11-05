import networkx as nx
from typing import Optional

class Arrow:
    def __init__(self, source : int, target : int, label : Optional[str]=None):
        self.source = source
        self.target = target
        self.label = label

    def __eq__(self, other : Arrow):
        return (self.source == other.source) and (self.target == other.target) and (self.label == other.label)

    def __repr__(self):
        return f"Arrow({self.source}->{self.target})"
    
class Path:
    def __init__(self,
                 arrows: tuple[Arrow,...] = (),
                 stationary_vertex: int = -1):
        
        if stationary_vertex != -1 and arrows:
            raise ValueError("Path cannot have both arrows and stationary_vertex.")

        self.stationary_vertex = stationary_vertex
        self.arrows = arrows

        # Check path is valid
        if self.arrows:
            for prev, curr in zip(self.arrows, self.arrows[1:]):
                if prev.target != curr.source:
                    raise ValueError(f"Invalid path: target of {prev} != source of {curr}")

    def is_stationary_path(self):
        return self.stationary_vertex is not None

    def source(self):
        if self.is_stationary_path():
            return self.stationary_vertex
        return self.arrows[0].source

    def target(self):
        if self.is_stationary_path():
            return self.stationary_vertex
        return self.arrows[-1].target
    
    def extend_at_end(self, arrow : Arrow):
        if self.target() == arrow.source:
            arrows = self.arrows + (arrow,)
            return Path(arrows)
    
    def extend_at_start(self, arrow : Arrow):
        if arrow.target == self.source():
            arrows = (arrow,) + self.arrows
            return Path(arrows)

    def __len__(self):
        if self.is_stationary_path():
            return 0
        return len(self.arrows)

    def __repr__(self):
        if self.is_stationary_path():
            return f"Stationary path at {self.stationary_vertex})"
        return f"Path = {self.arrows})"
    
class PathAlgebra:
    def __init__(self,
                 arrows : list[Arrow] = []):
        self.arrows = arrows

    def graph(self):
        return nx.MultiDiGraph([(arrow.source, arrow.target) for arrow in self.arrows])
    
    def vertices(self):
        return self.graph().nodes
    
    def stationary_paths(self):
        return [Arrow(vertex,vertex,"stationary") for vertex in self.vertices()]
    
    def is_path_of(self, path : Path):
        return all(path.arrows[i] in self.arrows for i in range(len(path)))
    
    def paths_from_fixed_pt_of_fixed_length(self, start : int, length : int):
        path = Path(stationary_vertex=start)
        paths = []
        def depth_first_search(path):
            if len(path) == length:
                paths.append(path)
                return
            for arrow in [a for a in self.arrows if a.source == path.target()]:
                new_path = path.extend_at_end(arrow)
                if new_path:
                    depth_first_search(new_path)
        depth_first_search(path)
        return paths
        
    def paths_from_fixed_pt(self, start : int, length : int = 8):
        paths = []
        for l in range(length+1):
            paths  += self.paths_from_fixed_pt_of_fixed_length( start, l)
        return paths
    
    def paths(self, length : int = 8):
        paths = {}
        for vertex in self.vertices():
            paths[vertex] = self.paths_from_fixed_pt(vertex, length)
        return paths

class MonomialQuiverAlgebra(PathAlgebra):
    def __init__(self,
                 arrows : list[Arrow] = [],
                 relations : list[Path] = [],
                 max_radical_length : int = 20):
        super().__init__(arrows)
        if not all(self.is_path_of(relations[i]) for i in range(len(relations))):
            raise ValueError(f"Invalid relations. {relations} must be a subset of {arrows}")
        self.relations = relations
        self.max_radical_length = max_radical_length

    def all_paths(self):
        return self.paths(self.max_radical_length)

    

    


if __name__ == "__main__":

    pa = PathAlgebra([Arrow(0,1),Arrow(1,2),Arrow(2,3),Arrow(3,0)])
    for vertex, paths in pa.paths().items():
        for path in paths:
            print(f"{vertex}: {path}")