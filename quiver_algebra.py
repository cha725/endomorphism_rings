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
                 stationary_vertex: int = None):
        
        if stationary_vertex is not None and arrows:
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
            if self.arrows is not None:
                arrows = self.arrows + (arrow,)
            else:
                arrows = (arrow,)
        return Path(arrows)
    
    def extend_at_start(self, arrow : Arrow):
        if arrow.target == self.source():
            if self.arrows is not None:
                arrows = (arrow,) + self.arrows
            else:
                arrows = (arrow,)
        return Path(arrows)
        
    def truncate(self, start_idx : int, end_idx : int):
        arrows = self.arrows[start_idx:end_idx]
        return Path(arrows)

    def is_subpath(self, other: Path):
        if len(self) == 0:
            return True
        if not(len(self) <= len(other)):
            return False
        if self.is_stationary_path():
            return self == Path(stationary_vertex=other.source())
        return any(other.truncate(i,i+len(self)) == self for i in range(len(other)-len(self)+1))

    def __len__(self):
        if self.is_stationary_path():
            return 0
        return len(self.arrows)
    
    def __eq__(self, other : Path):
        if self.is_stationary_path():
            return self.stationary_vertex == other.stationary_vertex
        return self.arrows == other.arrows

    def __repr__(self):
        if self.is_stationary_path():
            return f"Stationary path at {self.stationary_vertex}"
        return f"Path = {self.arrows}"
    
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

    def is_path(self, path : Path):
        return all(not relation.is_subpath(path) for relation in self.relations)

    def depth_first_search_paths(self, start: int, max_length: Optional[int] = None):
        if max_length is None:
            max_length = self.max_radical_length

        paths = [Path(stationary_vertex=start)]
        seen = []
        results = []
        while paths:
            path = paths.pop()
            if path in seen or not self.is_path(path):
                continue
            seen.append(path)
            results.append(path)
            if len(path) >= max_length:
                continue
            for arrow in [a for a in self.arrows if a.source == path.target()]:
                new_path = path.extend_at_end(arrow)
                if new_path:
                    paths.append(new_path)
        return results
    
    def paths(self, length : Optional[int] = None):
        paths = {}
        for vertex in self.vertices():
            paths[vertex] = self.depth_first_search_paths(vertex, length)
        return paths

    

    


if __name__ == "__main__":

    pa = PathAlgebra([Arrow(0,1),Arrow(1,2),Arrow(2,0)])

    qa = MonomialQuiverAlgebra(
    [Arrow(0,1),Arrow(1,2),Arrow(2,0)],
    [Path((Arrow(0,1),Arrow(1,2)))]
    )

    for vertex, paths in qa.paths().items():
        print(f"{vertex}: {paths}.")