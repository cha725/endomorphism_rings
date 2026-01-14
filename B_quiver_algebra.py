import networkx as nx

from A_bitmask_subgraph import Vertex, Arrow
    
class Path:
    """
    Represents a path in a directed graph.

    A path is either a stationary path at a single vertex or a tuple of composable
    arrows. Composability requires that for each consecutive pair of arrows,
    the target of the first matches the source of the second.

    Attributes:
        - arrows (tuple[Arrow] or Vertex):
            If a Vertex is provided, the path is the stationary path at that vertex.
            Otherwise, the path is the ordered concatenation of the list of arrows.
    """
    def __init__(self,
                 components: Vertex | tuple[Arrow,...]):
        
        self.vertex: Vertex | None = None
        self.arrows: tuple[Arrow,...] = tuple()

        if isinstance(components, Vertex):
            self.vertex = components
            self.arrows = tuple()
        else:
            self.vertex = None
            self.arrows = components

            if isinstance(self.arrows, tuple):
                for prev, curr in zip(self.arrows, self.arrows[1:]):
                    if prev.target != curr.source:
                        raise ValueError(f"Invalid path: target of {prev} != source of {curr}")

    def is_stationary_path(self) -> bool:
        return self.vertex is not None and not self.arrows

    def source(self):
        """ Returns source vertex of path. If stationary, then returns stationary vertex. """
        if not self.arrows:
            return self.vertex
        return self.arrows[0].source

    def target(self):
        """ Returns source vertex of path. If stationary, then returns stationary vertex. """
        if not self.arrows:
            return self.vertex
        return self.arrows[-1].target
    
    def vertices(self) -> list:
        """ Returns list of vertices in path. If stationary, then returns stationary vertex. """
        return [self.source()] + [a.target for a in self.arrows]
    
    def label(self) -> str:
        """ 
        Returns label of the path as concatenation of individual arrow labels. 
        For arrows a and b, the label ab means the path a then b.
        """
        if self.is_stationary_path():
            return f"Stationary path at {self.vertex}"
        return "".join(str(a.label) for a in self.arrows)
    
    def extend_at_end(self, arrow: Arrow) -> "Path | None":
        """ 
        If given arrow starts where path finishes, then returns concatenation path then arrow.
        Otherwise, returns None.
        """
        if self.target() != arrow.source:
            return None
        return Path(self.arrows + (arrow,))
    
    def extend_at_start(self, arrow: Arrow) -> "Path | None":
        """ 
        If given arrow ends where path starts, then returns new concatenation arrow then path.
        Otherwise, returns None.
        """
        if arrow.target != self.source():
            return None       
        return Path((arrow,) + self.arrows)
        
    def truncate(self, start_idx: int, end_idx: int) -> Path:
        """ Take subpath of path with indices [start_idx, ..., end_idx-1]. """
        arrows = self.arrows[start_idx:end_idx]
        return Path(arrows)

    def is_subpath(self, other: Path) -> bool:
        """ Returns True if self a subpath of other, otherwise False. """
        if len(self) == 0:
            return True
        if self.is_stationary_path():
            return self.source() == other.source()
        return any(other.truncate(i,i+len(self)) == self for i in range(len(other)-len(self)+1))
        
    def __len__(self):
        return len(self.arrows)
    
    def __eq__(self, other: Path):
        if not isinstance(other, Path):
            return False
        return (self.vertex == other.vertex) and (self.arrows == other.arrows)

    def __repr__(self):
        if self.is_stationary_path():
            return f"Path({self.label()})"
        vertices = "->".join(str(v.label) for v in self.vertices())
        return f"Path({self.label()}: {vertices})"
    
    def __hash__(self):
        if self.is_stationary_path():
            return hash(('stationary', self.vertex))
        return hash(tuple(self.arrows))

    
class PathAlgebra:
    """
    Path algebra of a directed graph.

    Attributes:
        - arrows (list[Arrow] or None):
            The list of directed edges of the graph. If None, the graph is
            treated as having no arrows.

        - vertices (list or None):
            If provided, these vertices are vertices of the graph. 
            The total set of vertices of the graph is this list plus any vertices
            appearing as sources or targets of the arrows.
            If None, the vertex set is the set of sources and targets of the arrows.

    Notes:
        - At least one of arrows or vertices must be non-empty; otherwise
            the path algebra would have no underlying graph.
        - The final vertex set is the union of the explicitly provided vertices
            (if any) and the vertices appearing in the arrows.
        - If the directed graph is not directed acyclic, then the path algebra
            will have infinitely many paths.
    """

    def __init__(self,
                 arrows: list[Arrow] | None = None,
                 vertices: list | None = None):

        self.arrows = arrows or []
        arrow_vertices: set[Vertex] = {a.source for a in self.arrows} | {a.target for a in self.arrows}        
        if vertices:
            arrow_vertices |= set(vertices)
        if not arrow_vertices:
            raise ValueError("Path algebra cannot be empty.")
        self.vertices: list[Vertex] = list(arrow_vertices)

    def graph(self) -> nx.MultiDiGraph:
        """ Return NetworkX graph representation of the path algebra """
        return nx.MultiDiGraph([(a.source, a.target) for a in self.arrows])
    
    def is_path_of(self, path: Path):
        """
        Check if the given path is a path in this path algebra.
        - Stationary paths are valid if their vertex is in the algebra.
        - Non-stationary paths are valid if all arrows are in the algebra.
        """
        if not path.arrows:
            return path.vertex in self.vertices
        return all(a in self.arrows for a in path.arrows)
    
# TODO: could make this take in gap code, then can use the quiver applet.

class MonomialQuiverAlgebra():
    def __init__(self,
                 arrows: Optional[list[Arrow]] = None,
                 relations: Optional[list[Path]] = None,
                 vertices: Optional[list] = None,
                 max_radical_length: int = 20):
        self.path_algebra = PathAlgebra(arrows, vertices)
        self.relations = relations or []
        if not all(self.is_path_of(r) for r in self.relations):
            raise ValueError(f"Invalid relations. {relations} must be a subset of {arrows}")
        self.max_radical_length = max_radical_length

    def is_path(self, path: Path):
        return all(not relation.is_subpath(path) for relation in self.relations)

    def is_path_extension_valid(self, path: Path, arrow: Arrow):
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
        initial_path = Path(vertex)
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
    
    def paths(self, max_length: Optional[int] = None, with_connections: bool = False) -> dict[int,list[Path]]:
        paths = {}
        for vertex in self.vertices:
            results, connections = self.dfs_paths_from_vertex(vertex, max_length)
            if with_connections:
                paths[vertex] = connections
            else:
                paths[vertex] = results
        return paths
    
    def __repr__(self):
        return f"MQA(vertices={self.vertices}, arrows={[a for a in self.arrows]}, relations={self.relations})"


###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################

if __name__ == "__main__":

    

    # TODO: Is this standard practice?

    class Examples:
        """
        Class to store examples of quivers.
        """
        def __init__(self,
                        examples: dict[str,MonomialQuiverAlgebra]):
            self.examples = examples

        def add(self, example: tuple[str, MonomialQuiverAlgebra]):
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
                        print(f"* {path}")

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
    
