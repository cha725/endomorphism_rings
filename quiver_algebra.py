import networkx as nx
from typing import Optional
from bitmask_subgraph import Arrow

    
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
    
    def vertices(self) -> list:
        if self.is_stationary_path():
            return [self.stationary_vertex]
        return [self.source()] + [a.target for a in self.arrows]
    
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
        if self.is_stationary_path():
            return self == Path(stationary_vertex=other.source())
        return any(other.truncate(i,i+len(self)) == self for i in range(len(other)-len(self)+1))
    
    def fancy_label(self):
        if self.is_stationary_path():
            return f"e({self.source})"
        return "".join(str(a.label) for a in self.arrows)
    
    def fancy_vertices(self):
        if self.is_stationary_path():
            return f"Stationary path at {self.source()}"
        return "->".join(str(v) for v in self.vertices())
        

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
        return f"{self.fancy_label()} : {self.fancy_vertices()}"
    
    def __hash__(self):
        if self.is_stationary_path():
            return hash(('stationary', self.stationary_vertex))
        return hash(tuple(self.arrows))

    
class PathAlgebra:
    def __init__(self,
                 arrows : Optional[list[Arrow]] = None,
                 vertices : Optional[list[int]] = None):
        
        if arrows is None and vertices is None:
            raise ValueError("Path algebra cannot be empty.")
        self.arrows = arrows or []
        self.arrow_vertices = {a.source for a in self.arrows} | {a.target for a in self.arrows}
        if vertices is None:
            self.vertices = list(self.arrow_vertices)
        else:
            self.vertices = list(self.arrow_vertices | set(vertices))

    def graph(self):
        return nx.MultiDiGraph([(a.source, a.target) for a in self.arrows])
    
    # TODO: dont like this - isn't mathematically accurate
    # def stationary_paths(self):
    #     return [Arrow(v,v,"stationary") for v in self.vertices]
    
    def is_path_of(self, path : Path):
        return all(a in self.arrows for a in path.arrows)
    
# TODO: could make this take in gap code, then can use the quiver applet.

class MonomialQuiverAlgebra(PathAlgebra):
    def __init__(self,
                 arrows : Optional[list[Arrow]] = None,
                 relations : Optional[list[Path]] = None,
                 vertices : Optional[list[int]] = None,
                 max_radical_length : int = 20):
        super().__init__(arrows, vertices)
        self.relations = relations or []
        if not all(self.is_path_of(r) for r in self.relations):
            raise ValueError(f"Invalid relations. {relations} must be a subset of {arrows}")
        self.max_radical_length = max_radical_length

    def is_path(self, path : Path):
        return all(not relation.is_subpath(path) for relation in self.relations)

    def is_path_extension_valid(self, path : Path, arrow : Arrow):
        """
        Given a valid path, check if the extension by an arrow is killed by the relations.
        """
        new_path = path.extend_at_end(arrow)
        if new_path is None:
            return False
        for rel in self.relations:
            if len(rel) <= len(new_path):
                if rel == new_path.truncate(len(new_path)-len(rel),len(new_path)):
                    return False
        return True


    def dfs_paths_from_vertex(self, vertex, max_length: Optional[int] = None) -> tuple[list[Path], list[Arrow]] | None:
        """
        Return list of paths of length at most max_length starting at given vertex.
        If max_length not given, then max_length set to be the max radical length of the quiver.
        """
        if vertex not in self.vertices:
            return None
        max_length = max_length or self.max_radical_length
        initial_path = Path(stationary_vertex=vertex)
        paths_to_check = [initial_path]
        connections: list[Arrow] = []
        seen = set()
        paths = []

        while paths_to_check:
            path = paths_to_check.pop()
            if path in seen:
                # if path seen move onto next path in paths
                continue
            # otherwise add to results and seen
            paths.append(path)
            seen.add(path)
            if len(path) >= max_length:
                # if path is maximum length move onto next path in paths
                continue
            # otherwise extend by arrows
            for arrow in [a for a in self.arrows if a.source == path.target()]:
                new_path = path.extend_at_end(arrow)
                if new_path and self.is_path(new_path):
                    paths_to_check.append(new_path)
                    connections.append(Arrow(path, 
                                             new_path, 
                                             arrow.label))
        return paths, connections
    
    def paths(self, length : Optional[int] = None, with_connections : bool = False) -> dict[int,list[Path]]:
        paths = {}
        for vertex in self.vertices:
            results, connections = self.dfs_paths_from_vertex(vertex, length)
            if with_connections:
                paths[vertex] = connections
            else:
                paths[vertex] = results
        return paths
    
    def __repr__(self):
        return f"MQA(vertices={self.vertices}, arrows={[a for a in self.arrows]}, relations={self.relations})"


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
            print(f"\n Vertices = {quiver.vertices}")
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
                  MonomialQuiverAlgebra(arrows=[Arrow(0,1,"a"),
                                                Arrow(1,2,"b"),
                                                Arrow(2,3,"c")],
                                        relations = [])))
    
    examples.add(("Type A rad2 relations",
                  MonomialQuiverAlgebra(arrows=[Arrow(0,1,"a"),
                                                Arrow(1,2,"b"),
                                                Arrow(2,3,"c")],
                                        relations = [Path((Arrow(0,1,"a"),
                                                           Arrow(1,2,"b"))),
                                                     Path((Arrow(1,2,"b"),
                                                           Arrow(2,3,"c")))])))

    examples.add(("Three cyclic",
                 MonomialQuiverAlgebra(arrows=[Arrow(0,1,"a"),
                                               Arrow(1,2,"b"),
                                               Arrow(2,0,"c")],
                                       relations=[Path((Arrow(0,1,"a"),
                                                        Arrow(1,2,"b")))])))

    examples.run()

    qa = MonomialQuiverAlgebra(arrows=[Arrow(0,1,"a"),
                                       Arrow(1,2,"b"),
                                       Arrow(2,3,"c")],
                               relations = [Path((Arrow(0,1,"a"),
                                                  Arrow(1,2,"b"))),
                                            Path((Arrow(1,2,"b"),
                                                  Arrow(2,3,"c")))])
    
