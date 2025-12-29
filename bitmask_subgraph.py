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
        self.vertex_to_index = {v: idx for idx, v in enumerate(self.vertices)}
        
        self.vertex_mask = (1 << self.num_vertices) - 1
        self.pred_mask = [0] * self.num_vertices
        self.succ_mask = [0] * self.num_vertices
        self.adj_mask = [0] * self.num_vertices
        for arrow in self.arrows:
            source_idx = self.vertex_to_index[arrow.source]
            target_idx = self.vertex_to_index[arrow.target]
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
        return [self._mask_to_vertices(mask) for mask in range(self.vertex_mask + 1)]

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
    def radical_layers(self) -> list[list[int]]:
        """
        Return radical layers as lists of vertices.
        """
        rad_layers = self._radical_layers
        return [self._mask_to_vertices(layer) for layer in rad_layers]

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
    def socle_layers(self) -> list[list[int]]:
        """
        Return socle layers as lists of vertices.
        """
        soc_layers = self._socle_layers
        return [self._mask_to_vertices(layer) for layer in soc_layers]

    def decompose(self, mask : int) -> list[list[int]]:
        """
        Return a list of connected components of the mask as lists of vertices.
        """
        if mask == 0:
            return [[]]
        # mask is connected iff starting at a single vertex can visit every other vertex
        
        summands = []
        remaining = mask
        while remaining:
            start = remaining & (~remaining + 1) # first non-zero bit in mask
            summand = start
            summand_remaining = start

            while summand_remaining:
            
                vertex = summand_remaining & (~summand_remaining + 1) # first non-zero bit
                summand_remaining &= ~vertex # remove vertex from remaining
                
                vertex_idx = vertex.bit_length()-1 # compute index of vertex
                
                neighbours = mask & self.adj_mask[vertex_idx] # adjacencies of vertex that are in the mask
                new_neighbours = neighbours & ~summand # find vertices that have not been visited
                
                summand |= new_neighbours # have visited the adjacent vertices now
                summand_remaining |= new_neighbours # add new neighbours to remaining

            summands.append(self._mask_to_vertices(summand))
            remaining &= ~summand # remove summand from remaining

        return summands
    
    def _is_connected(self, mask: int) -> bool:
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
    def compute_closed_subsets(self) -> tuple[list[list[int]],list[list[int]]]:
        """
        Compute which vertex subsets (represented as bitmasks) are closed under
        predecessors and/or successors.

        A subset S is:
            - predecessor-closed  iff  for every vertex v in S, all predecessors of v are in S.
            - successor-closed    iff  for every vertex v in S, all successors of v are in S.

        Only connected subsets are considered;
            disconnected masks are ignored and left out of the returned lists.

        Returns:
            tuple[list[list[int]],list[list[int]]]:
                - 0 entry = list of predecessor closed connected vertex subsets.
                - 1 entry = list of successor closed connected vertex subsets.
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

            pred_subsets = [self._mask_to_vertices(m) for m, val in enumerate(pred_results) if val]
            succ_subsets = [self._mask_to_vertices(m) for m, val in enumerate(succ_results) if val]

        return (pred_subsets, succ_subsets)    
    
    def compute_pred_closed_subsets(self) -> list[list[int]]:
        """
        Compute which vertex subsets (represented as bitmasks) are closed under
        predecessors.

        A subset S is predecessor-closed  iff  for every vertex v in S, all predecessors of v are in S.

        Only connected subsets are considered;
            disconnected masks are ignored and left out of the returned lists.

        Returns:
            list[list[int]]: list of predecessor closed connected vertex subsets.
        """        
        return self.compute_closed_subsets[0]
    
    def compute_succ_closed_subsets(self) -> list[list[int]]:
        """
        Compute which vertex subsets (represented as bitmasks) are closed under
        successors.

        A subset S is successor-closed    iff  for every vertex v in S, all successors of v are in S.

        Only connected subsets are considered;
            disconnected masks are ignored and left out of the returned lists.

        Returns:
            list[list[int]]: list of successor closed connected vertex subsets.
        """        
        return self.compute_closed_subsets[1]
    
    def radical_subgraphs(self):
        """
        Compute all successor-closed subgraphs generated by radical layers.
        """
        rad_layers = self._radical_layers
        subgraphs = []
        rad_subgraph = self.vertex_mask
        for layer in rad_layers:
            rad_subgraph &= ~layer
            ind_subgraphs = self.decompose(rad_subgraph)
            subgraphs.append(ind_subgraphs)
        return subgraphs
            



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
            rad_layers = graph.radical_layers
            soc_layers = graph.socle_layers
            layer_time = time.time() - start_time

            print(f"Time to compute radical and socle layers {layer_time:.4f} seconds")

            start_time = time.time()
            pred_closed = graph.compute_closed_subsets[0]
            succ_closed = graph.compute_closed_subsets[1]
            compute_time = time.time() - start_time

            print(f"Time to compute closed subsets: {compute_time:.4f} seconds")
            
            if verbose:

                print("\nVertices (with masks):")
                for vert, mask in enumerate(graph.vertex_to_index):
                    print(f" Vertex: {vert} -> Mask: {mask} = {bin(mask)} (binary)")
                
                print("\nArrows:")
                for arrow in graph.arrows:
                    print(f" {arrow}")
                
                print("\nRadical layers")
                for idx, layer in enumerate(rad_layers):
                    print(f" Radical layer: {idx} = Vertices: {layer}")

                print("\n Radical subgraphs")
                for idx, subgraphs in enumerate(graph.radical_subgraphs()):
                    print(f" Radical subgraph after removing layer {idx}:")
                    for sg in subgraphs:
                        print(f"  Vertices: {sg}")

                print("\nSocle layers")
                for idx, layer in enumerate(soc_layers):
                    print(f" Socle layer: {idx} = Vertices: {layer}")

                print("\nSubmodules:")
                for s in succ_closed:
                    print(f"Vertices: {s}")

                print("\nQuotients:")
                for p in pred_closed:
                    print(f"Vertices: {p}")

            times.append((init_time, layer_time, compute_time))

        return times

        

if __name__ == "__main__":

    examples = Examples({})

    for i in range(10):
        examples.add_random_graph(f"Random graph {i}", 5, 0.4)

    times = examples.run(verbose=True)

