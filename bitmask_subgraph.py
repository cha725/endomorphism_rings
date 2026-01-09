from typing import Optional
from functools import cached_property
import random, time

class Vertex:
    """
    Represents a vertex in a module diagram.

    Parameters:
        label: Label of the vertex.
        composition_factor: The label for the simple corresponding to this vertex.
        radical_layer (int): The value n such that the vertex is in the nth radical power
            of the module but not the (n+1)th radical power.
            Note such a value exists since the module diagram is finite.
            i.e. the maximum number of steps from a source to the vertex.
    """
    def __init__(self,
                 label,
                 composition_factor = None,
                 radical_layer: int | None = None):
        self.label = label
        self.composition_factor = composition_factor
        self.radical_layer = radical_layer

        if self.label is None:
            raise ValueError("Vertex label cannot be None.")

    def update_radical_layer(self, rad_layer: int):
        self.radical_layer = rad_layer

    def __eq__(self, other: "Vertex"):
        return self.label == other.label and self.composition_factor == other.composition_factor

    def __repr__(self):
        if self.composition_factor:
            return f"Vertex {self.label} with composition factor {self.composition_factor}"
        return f"Vertex {self.label}"
    
    def __hash__(self):
        return hash((self.label, self.composition_factor))

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
                 vertex_list: list[Vertex],
                 arrow_list : list[Arrow] | None):
        
        self.vertex_list = vertex_list
        self.num_vertices = len(self.vertex_list)

        self.vertex_to_index = {v : idx for idx, v in enumerate(self.vertex_list)}
        self.index_to_vertex = self.vertex_list

        self.vertex_mask = (1 << self.num_vertices) - 1
        self.vertex_to_mask = {v : (1 << idx) for v, idx in self.vertex_to_index.items()}
        self.mask_to_vertex = {(1 << idx) : v for v, idx in self.vertex_to_index.items()}

        self.arrow_list = arrow_list or []
        for arrow in self.arrow_list:
            if arrow.source not in self.vertex_list:
                raise ValueError(f"Source of {arrow} is {arrow.source} which is not a vertex {self.vertex_list}.")
            if arrow.target not in self.vertex_list:
                raise ValueError(f"Target of {arrow} is {arrow.target} which is not a vertex {self.vertex_list}.")
    
        self.num_arrows = len(self.arrow_list)

        self.pred_mask = [0] * self.num_vertices
        self.succ_mask = [0] * self.num_vertices
        self.adj_mask = [0] * self.num_vertices

        for arrow in self.arrow_list:
            source_idx = self.vertex_to_index[arrow.source]
            target_idx = self.vertex_to_index[arrow.target]
            self.succ_mask[source_idx] |= (1 << target_idx)
            self.pred_mask[target_idx] |= (1 << source_idx)
            self.adj_mask[source_idx] |= (1 << target_idx)
            self.adj_mask[target_idx] |= (1 << source_idx)

    ### VERTICES ###
    def _iterate_over_bits(self, mask: int) -> Iterator[int]:
        """ Returns interator that generates the individual bits in the given mask. """
        while mask:
            bit = mask & -mask
            mask ^= bit
            yield bit
    
    def _mask_to_indices(self, mask: int) -> list[int]:
        """ Converts mask to corresponding list of vertex indices. """
        remaining = mask
        indices = []
        while remaining:
            vertex = remaining & -remaining # get first non-zero bit
            remaining &= ~vertex # remove first non-zero bit
            vertex_idx = vertex.bit_length() - 1
            indices.append(vertex_idx)
        return indices
    
    def _indices_to_vertices(self, indices: list[int]) -> list[Vertex]:
        """ Converts list of vertex indices to corresponding list of vertices. """
        return [self.index_to_vertex[idx] for idx in indices]

    def _mask_to_vertices(self, mask: int) -> list[Vertex]:
        """ Converts mask to corresponding list of vertices. """
        return self._indices_to_vertices(self._mask_to_indices(mask))
    
    def _dim_of_mask(self, mask: int) -> int:
        """ Returns the number of vertices represented by the mask. """
        remaining = mask
        dim = 0
        while remaining:
            vertex = remaining & -remaining
            remaining &= ~vertex
            dim += 1
        return dim

    # Radical layers

    def sources_of_mask(self, mask: int) -> int:
        """
        Return bitmask of source vertices in the given mask.
        A source vertex is one with no predecessors in the mask.
        """
        sources = 0
        rem_mask = mask
        while rem_mask:
            v = rem_mask & -rem_mask
            rem_mask &= ~v
            idx = v.bit_length() - 1
            if self.pred_mask[idx] & mask == 0:
                sources |= v
        return sources
    
    @cached_property
    def _sources(self) -> int:
        """
        Return bitmask of vertices that are sources of the graph.
        A source vertex is one with no predecessors in the mask.
        """
        return self.sources_of_mask(self.vertex_mask)
    
    def sources(self) -> list[int]:
        """
        Return list of vertices that are sources.
        A source vertex is one with no predecessors in the mask.
        """
        return self._mask_to_vertices(self._sources)
    
    @cached_property
    def _radical_layers(self) -> list[int]:
        """
        Returns list of bitmasks representing radical layers.
        0th entry = sources, then removing sources gives next layer, etc.
        """
        remaining_mask = self.vertex_mask
        layers = []
        while remaining_mask:
            sources = self._sources_of_mask(remaining_mask)
            if sources == 0:
                return None
            layers.append(sources)
            remaining_mask &= ~sources
        return layers
    
    def radical_layers(self) -> list[list[Vertex]]:
        """
        Returns list of vertices representing radical layers.
        0th entry = sources, then removing sources gives next layer, etc.
        """
        layers = []
        if self._radical_layers:
            for layer_mask in self._radical_layers:
                layers.append(self._mask_to_vertices(layer_mask))
        return layers
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
    @cached_property
    def _connected_masks(self) -> list[bool]:
        """ 
        Returns list of booleans indexed by the masks. 
        The nth entry is true if the mask n is connected and False otherwise.
        """
        num_masks = 1 << self.num_vertices
        connected_masks = [False] * (num_masks)
        connected_masks[0] = True
        for mask in range(1, num_masks):
            first_bit = mask & -mask
            idx = first_bit.bit_length() - 1
            smaller_mask = mask ^ first_bit
            if smaller_mask == 0:
                # if smaller mask is 0 then mask is one vertex and hence connected.
                connected_masks[mask] = True
                continue
            if self.adj_mask[idx] & smaller_mask:
                # if first_bit is connected to a vertex in smaller_mask,
                # then mask is connected if and only if smaller_mask is.
                connected_masks[mask] = connected_masks[smaller_mask]
            else:
                # otherwise first_bit is not connected to the rest of the mask.
                connected_masks[mask] = False
        return connected_masks
    
    def _is_connected_mask(self, mask: int) -> bool:
        """
        Check if mask is connected (as an undirected graph).
        """
        return self._connected_masks[mask]
    
    def is_connected(self) -> bool:
        """ Check if the whole graph is connected. """
        return self._connected_masks[self.vertex_mask]

    def _connected_components_of(self, mask: int) -> list[int]:
        """ Returns a list of masks of connected components of the given mask. """
        if self._is_connected_mask(mask):
            return [mask]
        components = []
        remaining_mask = mask
        while remaining_mask:
            first_bit = remaining_mask & -remaining_mask
            component = first_bit
            to_check = first_bit
            while to_check:
                second_bit = to_check & -to_check
                to_check &= ~second_bit
                idx = second_bit.bit_length() - 1
                neighbours = mask & self.adj_mask[idx]
                new_neighbours = neighbours & ~component
                component |= new_neighbours
                to_check |= new_neighbours
            components.append(component)
            remaining_mask &= ~component
            if remaining_mask and self._is_connected_mask(remaining_mask):
                components.append(remaining_mask)
                break
        return components
    
    @cached_property
    def _connected_components(self) -> list[int]:
        """ Return list of connected components. """
        return self._connected_components_of(self.vertex_mask)
    
    @cached_property
    def connected_components(self) -> list[list[Vertex]]:
        """ Return list of connected components as lists of vertices. """
        return [self._mask_to_vertices(m) for m in self._connected_components]
    @cached_property
    def compute_closed_subsets(self) -> tuple[list[list[list[Vertex]]],list[list[list[Vertex]]]]:
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
        pred_subsets = [[] for _ in range(self.num_vertices + 1)]
        succ_subsets = [[] for _ in range(self.num_vertices + 1)]

        conn_comps = self._connected_components

        for comp in conn_comps:
            for mask in range(1, comp + 1):
                if not self._is_connected_mask(mask):
                    continue
                
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
                
                if pred_closed:
                    pred_results[mask] = True
                    dim = self._dim_of_mask(mask)
                    vertices = self._mask_to_vertices(mask)
                    if vertices not in pred_subsets[dim]:
                        pred_subsets[dim].append(vertices)
                if succ_closed:
                    succ_results[mask] = True
                    dim = self._dim_of_mask(mask)
                    vertices = self._mask_to_vertices(mask)
                    if vertices not in succ_subsets[dim]:
                        succ_subsets[dim].append(vertices)

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

