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

    def source(self) -> Vertex:
        """ Returns source vertex of path. If stationary, then returns stationary vertex. """
        if not self.arrows:
            return self.vertex
        return self.arrows[0].source

    def target(self) -> Vertex:
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
        
    def truncate(self, start_idx: int | None = None, end_idx: int | None = None) -> Path:
        """ Take subpath of path with indices [start_idx, ..., end_idx-1]. """
        new_arrows = self.arrows[start_idx:end_idx]
        if not new_arrows:
            return Path(self.source())
        return Path(new_arrows)
    
    def first_arrow(self) -> Arrow | None:
        """ Returns first arrow in the path or None is the path is stationary. """
        if self.is_stationary_path():
            return None
        return self.arrows[0]


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
        if self.vertex is None or other.vertex is None:
            if not (self.vertex is None and other.vertex is None):
                return False
        return (self.vertex == other.vertex) and (self.arrows == other.arrows)

    def __str__(self):
        if self.is_stationary_path():
            return f"e_{self.vertex}"
        return f"{self.label()}"

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
                 arrows: list[Arrow] | None = None,
                 relations: list[Path] | None = None,
                 vertices: list | None = None,
                 max_radical_length: int = 20):
        
        self.path_algebra = PathAlgebra(arrows, vertices)
        self.arrows = self.path_algebra.arrows
        self.vertices = self.path_algebra.vertices
        self.relations = relations or []
        if not all(self.path_algebra.is_path_of(r) for r in self.relations):
            raise ValueError(f"Invalid relations. {relations} must be a subset of {arrows}")
        self.max_radical_length = max_radical_length

    def is_path(self, path: Path) -> bool:
        """
        Check if the given path is a path in this quiver algebra.
        A path is valid if:
            - It is a path in the underlying path algebra.
            - It does not contain any of the monomial relations as a subpath.
        """
        if not self.path_algebra.is_path_of(path):
            return False
        return all(not relation.is_subpath(path) for relation in self.relations)

    def is_path_extension_valid(self, path: Path, arrow: Arrow) -> bool:
        """
        Given a path in the quiver algebra, check if the extension by the given arrow
        is also a path in the quiver algebra.
        
        Returns False if:
            - The concatenation of the path with the arrow is not a path.
            - The resulting concatenated path contains a relation as a subpath.
        """
        new_path = path.extend_at_end(arrow)
        if new_path is None:
            return False
        return self.is_path(new_path)


    def dfs_paths_from_vertex(self, 
                              vertex: Vertex, 
                              max_path_length: int | None = None
                              ) -> list[Path] | None:
        """
        Return list of paths of length at most max_path_length starting at given vertex.
        If max_path_length not given, then max_length set to be the max_radical_length of the quiver.
        Returns None if the vertex is not a vertex of the quiver algebra.
        """
        if vertex not in self.vertices:
            return None
        
        max_path_length = max_path_length or self.max_radical_length
        initial_path = Path(vertex)
        paths_to_check: list[Path] = [initial_path]
        seen: set[Path] = set()
        paths: list[Path] = []

        while paths_to_check:
            path = paths_to_check.pop()
            if path in seen:
                # if path seen move onto next path in paths
                continue
            # otherwise add to results and seen
            paths.append(path)
            seen.add(path)
            if len(path) >= max_path_length:
                # if path is maximum length move onto next path in paths
                continue
            # otherwise extend by arrows
            for arrow in [a for a in self.arrows if a.source == path.target()]:
                new_path = path.extend_at_end(arrow)
                if new_path and self.is_path(new_path):
                    paths_to_check.append(new_path)
        return paths
    
    def paths(self, max_path_length: int | None = None) -> dict[int,list[Path]]:
        """
        Returns dictionary with vertices as keys and a list of paths from that 
        vertex as values.
        """
        paths = {}
        max_path_length = max_path_length or self.max_radical_length
        for vertex in self.vertices:
            paths[vertex] = self.dfs_paths_from_vertex(vertex, max_path_length)
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
                  MonomialQuiverAlgebra(arrows=[Arrow(Vertex(0),Vertex(1),"a"),
                                                Arrow(Vertex(1),Vertex(2),"b"),
                                                Arrow(Vertex(2),Vertex(3),"c")],
                                        relations = [])))
    
    examples.add(("Type A rad2 relations",
                  MonomialQuiverAlgebra(arrows=[Arrow(Vertex(0),Vertex(1),"a"),
                                                Arrow(Vertex(1),Vertex(2),"b"),
                                                Arrow(Vertex(2),Vertex(3),"c")],
                                        relations = [Path((Arrow(Vertex(0),Vertex(1),"a"),
                                                           Arrow(Vertex(1),Vertex(2),"b"))),
                                                     Path((Arrow(Vertex(1),Vertex(2),"b"),
                                                           Arrow(Vertex(2),Vertex(3),"c")))])))

    examples.add(("Three cyclic",
                 MonomialQuiverAlgebra(arrows=[Arrow(Vertex(0),Vertex(1),"a"),
                                               Arrow(Vertex(1),Vertex(2),"b"),
                                               Arrow(Vertex(2),Vertex(0),"c")],
                                       relations=[Path((Arrow(Vertex(0),Vertex(1),"a"),
                                                        Arrow(Vertex(1),Vertex(2),"b")))])))

    examples.run()

    qa = MonomialQuiverAlgebra(arrows=[Arrow(Vertex(0),Vertex(1),"a"),
                                       Arrow(Vertex(1),Vertex(2),"b"),
                                       Arrow(Vertex(2),Vertex(3),"c")],
                               relations = [Path((Arrow(Vertex(0),Vertex(1),"a"),
                                                  Arrow(Vertex(1),Vertex(2),"b"))),
                                            Path((Arrow(Vertex(1),Vertex(2),"b"),
                                                  Arrow(Vertex(2),Vertex(3),"c")))])
    
