from typing import Optional
from functools import cached_property
import random, time

class Arrow:
    def __init__(self, source : int, target : int, label : Optional[str]=None):
        self.source = source
        self.target = target
        self.label = label

    def __eq__(self, other : Arrow):
        return (self.source == other.source) and (self.target == other.target) and (self.label == other.label)

    def __repr__(self):
        if self.label:
            return f"{self.label}: {self.source}->{self.target}"
        return f"{self.source}->{self.target}"
    
    def __hash__(self):
        return hash((self.source, self.target, self.label))

class BitmaskSubgraph:
    """
    Computes the predecessor and successor closed subsets of a directed graph.

    Uses bitmasks to compute these properties.
    For V vertices, the vertex bitmasks are on 0...V-1 bits.

    TODO: combine connectivity and closure functions for efficiency.
    """
    def __init__(self,
                 arrows : Optional[list[Arrow]] = None):
        if arrows is None:
            raise ValueError(f"Graph must contain at least one arrow.")
        self.arrows = arrows
        self.num_arrows = len(arrows)
        self.vertices = tuple(set(a.source for a in self.arrows) | set(a.target for a in self.arrows))
        self.index = {external_v: idx for idx, external_v in enumerate(self.vertices)}
        self.index_to_vertex = {v : k for k, v in self.index.items()}
        self.num_vertices = len(self.vertices)

    ### VERTICES ###

    @cached_property
    def vertex_mask(self) -> dict[int,int]:
        """
        Create bitmasks for each vertex and cache.
        For V vertices: vertex mask in bits 0..V-1.
        """
        return {external_v : 1 << idx for idx, external_v in enumerate(self.vertices)}

    @cached_property
    def full_vertex_mask(self) -> int:
        """
        Compute bitmask of all vertices.
        """
        return (1 << self.num_vertices) - 1
    
    @cached_property
    def mask_to_indices(self) -> list[list[int]]:
        """
        Create list to convert masks to the list of indices included.

        Returns:
            list: index i = list of vertex indices inside the binary expansion of i.
                    i.e. the indices of the 1s in the binary expansion.      
        """
        mask_to_vertex_indices = [[0]]
        for mask in range(1, self.full_vertex_mask + 1):
            remaining = mask
            subset = []
            while remaining:
                vertex = remaining & (~remaining + 1) # get first non-zero bit
                remaining &= ~vertex # remove first non-zero bit

                vertex_idx = vertex.bit_length() - 1

                subset.append(vertex_idx)
            mask_to_vertex_indices.append(subset)
        return mask_to_vertex_indices
    
    @cached_property
    def mask_to_vertices(self) -> list[list]:
        """
        Create list to convert masks to subsets of vertices.

        Returns:
            list: index i = list of vertices inside i.   
        """
        return [[self.vertices[idx] for idx in indices] for indices in self.mask_to_indices]
    
    @cached_property
    def adj_mask(self) -> list[dict[str,int]]:
        """
        Create predecessor "p" and successor "s" masks for each vertex.
        Returns:
            list: index i = {   "p" : predecessor vertex mask of vertex i, 
                                "s" : successor vertex mask of vertex i    }.
        """
        index = self.index
        adj = [{"p" : 0, "s" : 0} for _ in range(self.num_vertices)]
        for a in self.arrows:
            source_idx = index[a.source]
            target_idx = index[a.target]
            # a is an arrow a.source -> a.target
            # so a.source is a predecessor of a.target
            # and a.target is a successor of a.source
            adj[source_idx]["s"] |= (1 << target_idx)
            adj[target_idx]["p"] |= (1 << source_idx)
        return adj

    @cached_property
    def pred_mask(self) -> list[int]:
        """
        Create list of predecessor masks for each vertex.

        Returns:
            list: index i = predecessor vertex mask of vertex i.
        """
        return [adj_mask["p"] for adj_mask in self.adj_mask]

    @cached_property
    def succ_mask(self) -> list[int]:
        """
        Create list of successor masks for each vertex.

        Returns:
            list: index i = successor vertex mask of vertex i.
        """
        return [adj_mask["s"] for adj_mask in self.adj_mask]
    
    @cached_property
    def undirected_adj_mask(self) -> list[int]:
        """
        Create list of undirected neighbours for each vertex.

        Returns:
            list: index i = (undirected) neighbour vertex mask of vertex i.
        """
        return [adj_mask["p"] | adj_mask["s"] for adj_mask in self.adj_mask]
   
    def is_connected(self, mask : int) -> bool:
        """
        Check if mask is connected (as an undirected graph).
        """
        undir_adjs = self.undirected_adj_mask

        # mask is connected iff starting at a single vertex can visit every other vertex
        start = mask & (~mask + 1) # first non-zero bit in mask
        visited = start
        remaining = start

        while remaining:
            
            vertex = remaining & (~remaining + 1) # first non-zero bit
            remaining &= ~vertex # remove vertex from remaining
            
            vertex_idx = vertex.bit_length()-1 # compute index of vertex
            
            neighbours = mask & undir_adjs[vertex_idx] # adjacencies of vertex that are in the mask
            new_neighbours = neighbours & ~visited # find vertices that have not been visited
            
            visited |= new_neighbours # have visited the adjacent vertices now
            remaining |= new_neighbours # add new neighbours to remaining

        return visited == mask
    
    @cached_property
    def closed_subsets(self) -> tuple[list[bool],...]:
        """
        Compute which vertex subsets (represented as bitmasks) are closed under
        predecessors and/or successors.

        A subset S is:
            - predecessor-closed  iff  for every vertex v in S, all predecessors of v are in S.
            - successor-closed    iff  for every vertex v in S, all successors of v are in S.

        Only connected subsets are considered;
            disconnected masks are ignored and left out of the returned lists.

        Returns:
            dict:   {   "pred_closed" : list[bool],
                            A boolean list indexed by mask i for i in (0,full_mask].
                            Entry i is True iff mask i is connected and predecessor-closed.

                        "succ_closed" : list[bool],
                            A boolean list indexed by mask i for i in (0,full_mask].
                            Entry i is True iff mask i is connected and successor-closed.
                    }   

        Notes:
            - the 0 index entry of the list is 0 (both closed under predecessors and successors)
            - if the mask is disconnected, its entry in both lists is False
            - masks correspond to subsets of vertices
    """        
        full_mask = self.full_vertex_mask
        mask_to_indices = self.mask_to_indices
        pred_mask = self.pred_mask
        succ_mask = self.succ_mask
        pred_closed = [False for _ in range(full_mask + 1)]
        succ_closed = [False for _ in range(full_mask + 1)]

        for mask in range(1, full_mask + 1):
            if not self.is_connected(mask):
                continue

            p_closed = True
            s_closed = True

            for vertex_idx in mask_to_indices[mask]:
            
                p_mask = pred_mask[vertex_idx]
                s_mask = succ_mask[vertex_idx]
                    
                if p_closed and p_mask & (~mask): # check if p_mask outside of mask
                    p_closed = False
                if s_closed and s_mask & (~mask): # check if s_mask outside of mask
                    s_closed = False
                
                if not(p_closed or s_closed):
                    break

            pred_closed[mask] = p_closed
            succ_closed[mask] = s_closed

        return (pred_closed, succ_closed)
    
    @cached_property
    def pred_closed_subsets(self):
        return [idx for idx, is_closed in enumerate(self.closed_subsets[0]) if is_closed]
    
    @cached_property
    def succ_closed_subsets(self):
        return [idx for idx, is_closed in enumerate(self.closed_subsets[1]) if is_closed]

class Examples:
    """
    Class to store examples of bitmask graphs.
    """
    def __init__(self, 
                 examples: dict[str, list[Arrow]]):
        self.examples = examples

    def add(self, name: str, arrows: list[Arrow]):
        self.examples[name] = arrows

    def add_random_graph(self, name: str, num_vertices: int, edge_prob: float):
        """
        Add a random directed tree with num_vertices vertices and edge probability edge_prob.
        """
        arrows = []
        for i in range(num_vertices):
            for j in range(i+1,num_vertices):
                if random.random() < edge_prob:
                    arrows.append(Arrow(i, j))
        self.add(name, arrows)

    def run(self):
        times = []
        for name, arrows in self.examples.items():
            print(f"\n=== BitmaskGraph Example: {name} ===")
            
            start_time = time.time()
            graph = BitmaskSubgraph(arrows)
            init_time = time.time() - start_time

            print(f"\nTime to initialise graph: {init_time:.4f} seconds")

            start_time = time.time()
            pred_closed = graph.pred_closed_subsets
            succ_closed = graph.succ_closed_subsets
            compute_time = time.time() - start_time

            print(f"Time to compute closed subsets: {compute_time:.4f} seconds")
            
            print("\nVertices (with masks):")
            for vert, mask in graph.vertex_mask.items():
                print(f" Vertex: {vert} -> Mask: {mask} = {bin(mask)} (binary)")

            mask_to_vertices = graph.mask_to_vertices

            print("\nSubmodules:")
            for s in succ_closed:
                print(f" Mask: {s} = {bin(s)} -> Vertices: {mask_to_vertices[s]}")

            print("\nQuotients:")
            for q in pred_closed:
                print(f" Mask: {q} = {bin(q)} -> Vertices: {mask_to_vertices[q]}")

            times.append((init_time, compute_time))

        return times

        

if __name__ == "__main__":

    examples = Examples({})

    # examples.add("Type A", [Arrow(0,1), Arrow(1,2), Arrow(2,3)])

    # examples.add("Diamond", [Arrow(0,1), Arrow(0,2), Arrow(1,3), Arrow(2,3)])

    for i in range(10):
        examples.add_random_graph(f"Random graph {i}", 20, 0.25)

    times = examples.run()
    for idx, t in enumerate(times):
        print(f"\n=== Example {idx} ===")
        print(f"\n Initialisation time {t[0]:.4f} seconds; \n Computation time {t[1]:.4f} seconds.")



