import networkx as nx
from typing import Optional
from bitmaskgraph import Arrow

    
class Path:
    def __init__(self,
                 arrows: tuple[Arrow,...] = (),
                 stationary_vertex: Optional[int] = None):
        
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
        if self.target() != arrow.source:
            return None
        return Path(self.arrows + (arrow,))
    
    def extend_at_start(self, arrow : Arrow):
        if arrow.target != self.source():
            return None
        return Path((arrow,) + self.arrows)
        
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
        return f"Path{self.arrows}"
    
    def __hash__(self):
        if self.is_stationary_path():
            return hash(('stationary', self.stationary_vertex))
        return hash(tuple(self.arrows))

    
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
        connecting_edges = [(Path(stationary_vertex=start), None, Path(stationary_vertex=start))]
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
                    connecting_edges.append((path, arrow, new_path))
        return results, connecting_edges
    
    def paths(self, length : Optional[int] = None, with_connections : bool = False) -> dict[int,list]:
        paths = {}
        for vertex in self.vertices():
            if with_connections:
                paths[vertex] = self.depth_first_search_paths(vertex, length)[1]
            else:
                paths[vertex] = self.depth_first_search_paths(vertex, length)[0]
        return paths

    


# TODO: Is this standard practice?

class Examples:
    """
    Class to store examples of quivers.
    """
    def __init__(self,
                    examples : dict[str,MonomialQuiverAlgebra]):
        self.examples = examples

    def add(self, example : tuple[str, MonomialQuiverAlgebra]):
        self.examples[example[0]] = example[1]

    def run(self):
        for name, quiver in self.examples.items():
            print(f"\n=== Example: {name} ===")
            print(f"\n Vertices = {quiver.vertices()}")
            print(f" Arrows = {quiver.arrows}")
            print(f" Relations = {quiver.relations}")
            print(f"\n-- Paths in quiver --")
            for vertex, paths in quiver.paths().items():
                print(f"\n- Starting at vertex {vertex} -")
                for path in paths:
                    print(f"• {path}")
                      


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

    examples.run()