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
                 arrows : Optional[tuple[Arrow,...]] = None):
        if arrows is None:
            raise ValueError(f"Graph must contain at least one arrow.")
        self.arrows = arrows
        self.num_arrows = len(arrows)
        self.vertices = tuple(sorted(set(a.source for a in self.arrows) | set(a.target for a in self.arrows)))
        self.num_vertices = len(self.vertices)
        self.index = {v: idx for idx, v in enumerate(self.vertices)}
        self.index_to_vertex = {v : k for k, v in self.index.items()}
        
        self.vertex_mask = (1 << self.num_vertices) - 1
        self.pred_mask = [0] * self.num_vertices
        self.succ_mask = [0] * self.num_vertices
        self.adj_mask = [0] * self.num_vertices
        for arrow in self.arrows:
            source_idx = self.index[arrow.source]
            target_idx = self.index[arrow.target]
            self.succ_mask[source_idx] |= (1 << target_idx)
            self.pred_mask[target_idx] |= (1 << source_idx)
            self.adj_mask[source_idx] |= (1 << target_idx)
            self.adj_mask[target_idx] |= (1 << source_idx)
        # TODO: does pred_mask etc. suggest an integer rather than a list?

    ### VERTICES ###
    
    def _mask_to_vertices(self, mask: int) -> list:
        """
        Returns list of vertices in given mask.
        """
        remaining = mask
        vertices = []
        while remaining:
            vertex = remaining & (~remaining + 1) # get first non-zero bit
            remaining &= ~vertex # remove first non-zero bit
            vertex_idx = vertex.bit_length() - 1
            vertices.append(vertex_idx)
        return vertices

    @cached_property
    def mask_to_vertices(self) -> list[list]:
        """
        Create list to convert masks to subsets of vertices.

        Returns:
            list: index i = list of vertices inside i.   
        """
        return [self._mask_to_vertices(mask) for mask in range(self.full_vertex_mask + 1)]

    @cached_property
    def sources(self) -> list[bool]:
        """
        Create list of vertices that are sources.
            i.e. No predecessors.

        Returns:
            list: index i = True if vertex i is a source.
        """
        return [pred == 0 for pred in self.pred_mask]

    @cached_property
    def sinks(self) -> list[bool]:
        """
        Create list of vertices that are sinks.
            i.e. No successors.

        Returns:
            list: index i = True if vertex i is a sink.
        """
        return [succ == 0 for succ in self.succ_mask]




    @cached_property
    def _radical_layers(self) -> list[int]:
        """
        Compute radical layers.
        A vertex has radical layer n if the longest path from any source to it has length n.

        Returns:
            list: index i = vertices in radical layer i.
        """
        succ_mask = self.succ_mask
        radical_layers = []
        visited = 0
        current_layer = 0

        for idx, is_source in enumerate(self.sources):
            if is_source:
                current_layer |= 1 << idx

        while current_layer:
            radical_layers.append(current_layer)
            visited |= current_layer

            remaining = current_layer
            new_layer = 0

            while remaining:
                v = remaining & -remaining
                remaining &= ~v
                v_idx = v.bit_length() - 1
                succ = succ_mask[v_idx]
                while succ:
                    succ_v = succ & -succ
                    succ &= ~succ_v
                    succ_idx = succ_v.bit_length() - 1
                    if self.pred_mask[succ_idx] & ~visited == 0:
                        new_layer |= succ_v
                
            current_layer = new_layer & ~visited # remove already visited vertices

        return radical_layers
    
    @cached_property
    def radical_layer_vertices(self):
        """
        Create list of radical layers as lists of vertices.
            i.e. minimum number of steps from a source.

        Returns:
            list: index i = all vertices in radical layer i.
        """
        return [self.mask_to_vertices[layer] for layer in self.radical_layer_mask]

    @cached_property
    def _socle_layers(self) -> list[int]:
        """
        Create list of socle layers as bitmasks.
            i.e. maximum number of steps from a sink.

        Returns:
            list: index i = vertices in socle layer i.
        """
        socle_layers = []
        pred_mask = self.pred_mask
        next_layer = 0
        for idx, is_sink in enumerate(self.sinks):
            if is_sink:
                next_layer |= 1 << idx
            
        visited = next_layer
        while next_layer:
            socle_layers.append(next_layer)

            remaining = next_layer
            new_layer = 0
            while remaining:
                v = remaining & -remaining
                remaining &= ~v
                v_idx = v.bit_length() - 1
                pred = pred_mask[v_idx]
                while pred:
                    pred_v = pred & -pred
                    pred &= ~pred_v
                    pred_idx = pred_v.bit_length() - 1
                    if self.succ_mask[pred_idx] & ~visited == 0:
                        new_layer |= pred_v
 
            next_layer = new_layer &~ visited
            visited |= next_layer

        return socle_layers
    
    @cached_property
    def socle_layer_vertices(self):
        """
        Create list of socle layers.

        Returns:
            list: index i = all vertices in socle layer i.
        """
        return [self.mask_to_vertices[layer] for layer in self.socle_layers_mask]

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
        if mask == 0:
            return True
        # mask is connected iff starting at a single vertex can visit every other vertex
        start = mask & (~mask + 1) # first non-zero bit in mask
        visited = start
        remaining = start

        while remaining:
            
            vertex = remaining & (~remaining + 1) # first non-zero bit
            remaining &= ~vertex # remove vertex from remaining
            
            vertex_idx = vertex.bit_length()-1 # compute index of vertex
            
            neighbours = mask & self.adj_mask[vertex_idx] # adjacencies of vertex that are in the mask
            new_neighbours = neighbours & ~visited # find vertices that have not been visited
            
            visited |= new_neighbours # have visited the adjacent vertices now
            remaining |= new_neighbours # add new neighbours to remaining

        return visited == mask
    
    @cached_property
    def compute_closed_subsets(self) -> tuple[list[bool],list[bool]]:
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
        full_mask = self.vertex_mask
        pred_mask = self.pred_mask
        succ_mask = self.succ_mask

        pred_results = [False for _ in range(full_mask + 1)]
        succ_results = [False for _ in range(full_mask + 1)]

        for mask in range(1, full_mask + 1):
            connected = False
            pred_closed = True
            succ_closed = True

            for vertex_idx in range(self.num_vertices):
                if not (mask & (1 << vertex_idx)):
                    continue # vertex not in mask
                if pred_mask[vertex_idx] & (~mask): # check if p_mask outside of mask
                    pred_closed = False
                if succ_mask[vertex_idx] & (~mask): # check if s_mask outside of mask
                    succ_closed = False
                if not(pred_closed or succ_closed):
                    break # no need to continue if neither closed

            if pred_closed or succ_closed:
                if self._is_connected(mask):
                    connected = True

            pred_results[mask] = pred_closed and connected
            succ_results[mask] = succ_closed and connected

        return (pred_results, succ_results)
    
    @cached_property
    def connected_closed_subsets(self):
        """
        Compute connected predecessor-closed and successor-closed subsets efficiently.

    #     TODO: Why is this slower than the original function?
    #     """
    #     n = len(self.vertex_mask)
    #     full_mask = self.full_vertex_mask
    #     pred = self.pred_mask
    #     succ = self.succ_mask
    #     undir = self.undirected_adj_mask

    #     pred_missing = [0] * (full_mask + 1)
    #     succ_missing = [0] * (full_mask + 1)
    #     pred_closed = [False] * (full_mask + 1)
    #     succ_closed = [False] * (full_mask + 1)
    #     connected_comp = [0] * (full_mask + 1)

    #     pred_results = [False] * (full_mask + 1)
    #     succ_results = [False] * (full_mask + 1)

    #     # mask 0 is trivially closed and connected
    #     pred_closed[0] = succ_closed[0] = True
    #     connected_comp[0] = 0
    #     pred_results[0] = succ_results[0] = True

    #     # Process masks in increasing number of vertices in mask
    #     masks_by_size = [[] for _ in range(n + 1)]
    #     for mask in range(1, full_mask + 1):
    #         size = bin(mask).count("1")
    #         masks_by_size[size].append(mask)

    #     for size in range(1, n + 1):
    #         for mask in masks_by_size[size]:
    #             v = mask & -mask
    #             prev_mask = mask & (mask - 1)
    #             v_idx = v.bit_length() - 1

    #             # Check predecessor / successor closure
    #             p_missing = (pred[v_idx] & ~mask) | (pred_missing[prev_mask] & ~v)
    #             s_missing = (succ[v_idx] & ~mask) | (succ_missing[prev_mask] & ~v)
    #             pred_missing[mask] = p_missing
    #             succ_missing[mask] = s_missing

    #             if p_missing == 0 and (pred_closed[prev_mask] or pred_missing[prev_mask] == v):
    #                 pred_closed[mask] = True
    #             if s_missing == 0 and (succ_closed[prev_mask] or succ_missing[prev_mask] == v):
    #                 succ_closed[mask] = True
    #             if not pred_closed[mask] and not succ_closed[mask]:
    #                 continue  # move on if neither pred or succ closed

    #             # Connectivity
    #             prev_comp = connected_comp[prev_mask]
    #             if prev_mask == 0:
    #                 new_comp = mask
    #             else:
    #                 # Only check connectivity if mask is potentially closed
    #                 if prev_comp & undir[v_idx]:
    #                     # BFS to expand component
    #                     visited = prev_comp | v
    #                     remaining = visited
    #                     while remaining:
    #                         bit = remaining & -remaining
    #                         remaining ^= bit
    #                         idx = bit.bit_length() - 1
    #                         neighbors = (mask & undir[idx]) & ~visited
    #                         visited |= neighbors
    #                         remaining |= neighbors
    #                     new_comp = visited
    #                 else:
    #                     new_comp = prev_comp  # disconnected

    #             connected_comp[mask] = new_comp
    #             is_conn = new_comp == mask

    #             pred_results[mask] = pred_closed[mask] and is_conn
    #             succ_results[mask] = succ_closed[mask] and is_conn

    #     return pred_results, succ_results
    

    @cached_property
    def pred_closed_index_subsets(self):
        return [idx for idx, is_closed in enumerate(self.closed_subsets[0]) if is_closed]
    
    @cached_property
    def succ_closed_index_subsets(self):
        return [idx for idx, is_closed in enumerate(self.closed_subsets[1]) if is_closed]
    
    @cached_property
    def pred_closed_vertex_subsets(self):
        return [self.mask_to_vertices[mask] for mask in self.pred_closed_index_subsets]
    
    @cached_property
    def succ_closed_vertex_subsets(self):
        return [self.mask_to_vertices[mask] for mask in self.succ_closed_index_subsets]
    

class Examples:
    """
    Class to store examples of bitmask graphs.
    """
    def __init__(self, 
                 examples: dict[str, tuple[Arrow]]):
        self.examples = examples

    def add(self, name: str, arrows: tuple[Arrow]):
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
        self.add(name, tuple(arrows))

    def run(self, verbose : bool = False):
        times = []
        for name, arrows in self.examples.items():
            print(f"\n=== BitmaskGraph Example: {name} ===")
            
            start_time = time.time()
            graph = BitmaskSubgraph(arrows)
            init_time = time.time() - start_time

            print(f"\nTime to initialise graph: {init_time:.4f} seconds")

            start_time = time.time()
            rad_layers = graph.radical_layer_vertices
            soc_layers = graph.socle_layer_vertices
            layer_time = time.time() - start_time

            print(f"Time to compute radical or socle layers {layer_time:.4f} seconds")

            start_time = time.time()
            pred_closed = graph.pred_closed_index_subsets
            succ_closed = graph.succ_closed_index_subsets
            compute_time = time.time() - start_time

            print(f"Time to compute closed subsets: {compute_time:.4f} seconds")
            
            if verbose:

                print("\nVertices (with masks):")
                for vert, mask in graph.vertex_mask.items():
                    print(f" Vertex: {vert} -> Mask: {mask} = {bin(mask)} (binary)")
                
                print("\nRadical layers")
                for idx, layer in enumerate(rad_layers):
                    print(f" Radical layer: {idx} = Vertices: {layer}")

                print("\nSocle layers")
                for idx, layer in enumerate(soc_layers):
                    print(f" Socle layer: {idx} = Vertices: {layer}")

                print("\nSubmodules:")
                for s in succ_closed:
                    print(f" Mask: {s} = {bin(s)} -> Vertices: {graph.mask_to_indices[s]}")

                print("\nQuotients:")
                for q in pred_closed:
                    print(f" Mask: {q} = {bin(q)} -> Vertices: {graph.mask_to_indices[s]}")

            times.append((init_time, layer_time, compute_time))

        return times

        

if __name__ == "__main__":

    examples = Examples({})

    for i in range(10):
        examples.add_random_graph(f"Random graph {i}", 10, 0.4)

    times = examples.run(verbose=True)

