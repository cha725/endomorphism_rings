from functools import cached_property
from typing import Iterator
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
    """
    Represents an arrow as a pair (source, target) with an optional label.
    """
    def __init__(self, 
                 source, 
                 target, 
                 label = None):
        self.source = source
        self.target = target
        self.label = label

    def __eq__(self, other: "Arrow") -> bool:
        return (self.source == other.source) and (self.target == other.target) and (self.label == other.label)

    def __repr__(self):
        if self.label:
            return f"Arrow({self.label}: {self.source} -> {self.target})"
        return f"Arrow({self.source} -> {self.target})"
    
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
    
    def sources(self) -> list[Vertex]:
        """
        Return list of vertices that are sources.
        A source vertex is one with no predecessors in the mask.
        """
        return self._mask_to_vertices(self._sources)
    
    @cached_property
    def _radical_layers(self) -> list[int] | None:
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
    
    def _compute_radical_subgraphs(self) -> list[int]:
        """ Returns list of bitmasks where ith entry is all the vertices in radical layer j>=i. """
        rad_layers = self._radical_layers
        subgraphs = []
        if rad_layers:
            sub_rad_layers = rad_layers.copy()
            for layer_idx, layer in enumerate(rad_layers):
                subgraphs.append(layer)
                sub_rad_layers.remove(layer)
                for sub_layers in sub_rad_layers:
                    subgraphs[layer_idx] |= sub_layers
        return subgraphs
    
    def compute_radical_subgraphs(self) -> list[list[Vertex]]:
        """ Returns list of bitmasks where ith entry is all the vertices in radical layer j>=i. """
        subgraphs = self._compute_radical_subgraphs()
        ind_subgraphs = []
        for subgraph in subgraphs:
            summands = [self._mask_to_vertices(m) for m in self._connected_components_of(subgraph)]
            ind_subgraphs += summands
        return ind_subgraphs

    # Socle layers

    def _sinks_of_mask(self, mask: int) -> int:
        """
        Return bitmask of sink vertices in the given mask.
        A sink vertex is one with no successors in the mask.
        """
        sinks = 0
        for bit in self._iterate_over_bits(mask):
            idx = bit.bit_length() - 1
            if self.succ_mask[idx] & mask == 0:
                sinks |= bit
        return sinks
    
    @cached_property
    def _sinks(self) -> int:
        """
        Return bitmask of vertices that are sinks of the graph.
        A sink vertex is one with no successors in the mask.
        """
        return self._sinks_of_mask(self.vertex_mask)
    
    def sinks(self) -> list[Vertex]:
        """
        Return list of sink vertices.
        """
        return self._mask_to_vertices(self._sinks)
    
    @cached_property
    def _socle_layers(self) -> list[int]:
        """
        Returns list of bitmasks representing socle layers.
        0th entry = sinks, then removing sinks gives next layer, etc.
        """
        remaining_mask = self.vertex_mask
        layers = []
        while remaining_mask:
            sinks = self._sinks_of_mask(remaining_mask)
            layers.append(sinks)
            remaining_mask &= ~sinks
        return layers
    
    def socle_layers(self) -> list[list[Vertex]]:
        """
        Returns list of vertices representing socle layers.
        0th entry = sinks, then removing sinks gives next layer, etc.
        """
        layers = []
        for layer_mask in self._socle_layers:
            layers.append(self._mask_to_vertices(layer_mask))
        return layers
    
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
    
    def _is_pred_closed(self, mask: int) -> bool:
        """
        Check if a bitmask is closed under predecessors.
        Returns True if closed under predecessors, or False otherwise.
        TODO: would it be faster to compute the sources first?
        """
        rem_mask = mask
        while rem_mask:
            v = rem_mask & -rem_mask
            rem_mask &= ~v
            v_idx = v.bit_length() - 1
            pred_v = self.pred_mask[v_idx]
            rem_pred_v = pred_v
            while rem_pred_v:
                u = rem_pred_v & -rem_pred_v
                rem_pred_v &= ~u
                if u & mask == 0:
                    return False
        return True

    
    def _is_succ_closed(self, mask: int) -> bool:
        """
        Check if a bitmask is closed under successors.
        Returns True if closed under successors, or False otherwise.
        TODO: would it be faster to compute the sinks first?
        """
        rem_mask = mask
        while rem_mask:
            v = rem_mask & -rem_mask
            rem_mask &= ~v
            v_idx = v.bit_length() - 1
            succ_v = self.succ_mask[v_idx]
            rem_succ_v = succ_v
            while rem_succ_v:
                u = rem_succ_v & -rem_succ_v
                rem_succ_v &= ~u
                if u & mask == 0:
                    return False
        return True
      

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
    
    def compute_pred_closed_subsets(self) -> list[list[list[Vertex]]]:
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
    
    def compute_succ_closed_subsets(self) -> list[list[list[Vertex]]]:
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
            vertices = [Vertex(str(i)) for i in range(num_vertices)]
            label_to_vertex = {v.label : v for v in vertices}
            arrows = []
            for i in range(num_vertices):
                for j in range(i+1,num_vertices):
                    if random.random() < edge_prob:
                        arrows.append(Arrow(label_to_vertex[str(i)], label_to_vertex[str(j)]))
            self.add(name, arrows)

        def run(self, verbose : bool = False):
            times = []
            for name, arrows in self.examples.items():
                print(f"\n=== BitmaskGraph Example: {name} ===")
                
                vertices = list(set(a.source for a in arrows) | set(a.target for a in arrows))

                start_time = time.time()
                graph = BitmaskSubgraph(vertices,arrows)
                init_time = time.time() - start_time


                print(f"\nTime to initialise graph: {init_time:.4f} seconds")

                start_time = time.time()
                rad_layers = graph.radical_layers()
                soc_layers = graph.socle_layers()
                layer_time = time.time() - start_time

                print(f"Time to compute radical and socle layers {layer_time:.4f} seconds")

                start_time = time.time()
                rad_subgraphs = graph.compute_radical_subgraphs()
                rad_sub_time = time.time() - start_time

                print(f"Time to compute radical subgraphs: {rad_sub_time:.4f} seconds")

                start_time = time.time()
                pred_closed = graph.compute_closed_subsets[0]
                succ_closed = graph.compute_closed_subsets[1]
                compute_time = time.time() - start_time

                print(f"Time to compute closed subsets: {compute_time:.4f} seconds")
                

                if verbose:

                    print("\nVertices (with masks):")
                    for vert, mask in graph.vertex_to_index.items():
                        print(f" {vert} -> Mask: {mask}")
                    
                    print("\nArrows:")
                    for arrow in graph.arrow_list:
                        print(f" {arrow}")
                    
                    print("\nRadical layers")
                    for idx, layer in enumerate(rad_layers):
                        print(f" Radical layer: {idx} = Vertices: {[v.label for v in layer]}")

                    print("\nRadical graphs")
                    for rad_graph in rad_subgraphs:
                        print(f"Vertices: {[v.label for v in rad_graph]}")

                    print("\nSocle layers")
                    for idx, layer in enumerate(soc_layers):
                        print(f" Socle layer: {idx} = Vertices: {[v.label for v in layer]}")

                    print("\nConnected components:")
                    for c in graph.connected_components:
                        print(f"Vertices: {[v.label for v in c]}")

                    print("\nSubmodules:")
                    for list_s in succ_closed:
                        for s in list_s:
                            print(f"Vertices: {[v.label for v in s]}")

                    print("\nQuotients:")
                    for list_p in pred_closed:
                        for p in list_p:
                            print(f"Vertices: {[v.label for v in p]}")

                times.append((init_time, layer_time, compute_time))

            return times

            


    examples = Examples({})

    for i in range(1):
        examples.add_random_graph(f"Random graph {i}", 15, 0.4)

    times = examples.run(verbose=True)
