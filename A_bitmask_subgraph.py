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

    Attributes:
        - source: The start of the arrow.
        - target: The end of the arrow.
        - label: The label of the arrow.
    """
    def __init__(self, 
                 source: Vertex, 
                 target: Vertex, 
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

class BitmaskGraph:
    """
    Represents a directed acyclic graph with vertices encoded as bitmasks.

    Each vertex is assigned a bitmask corresponding to its index in `vertex_list`.
    These bitmasks allow efficient computation of predecessor masks, successor
    masks, adjacency, radical layers, and connectivity/closure properties.

    The graph is assumed to be a directed acyclic graph (DAG) and thus has no
    directed cycles, but may have undirected cycles.

    Attributes:
        vertex_list (list[Vertex]):
            The list of vertices in the graph.
        arrow_list (list[Arrow]):
            The directed edges (arrows) of the graph. 
            The source and target of each arrow must be a member of vertex_list.
        pred_masks (list[int]):
            For each vertex index i, a bitmask of all predecessors of i.
        succ_masks (list[int]):
            For each vertex index i, a bitmask of all successors of i.
        adj_masks (list[int]):
            For each vertex index i, a bitmask of all adjacent vertices
            (treating edges as undirected).

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
        self.num_arrows = len(self.arrow_list)

        self.pred_masks = [0] * self.num_vertices
        self.succ_masks = [0] * self.num_vertices
        self.adj_masks = [0] * self.num_vertices

        for arrow in self.arrow_list:
            try:
                source_idx = self.vertex_to_index[arrow.source]
            except:
                raise ValueError(f"Source of {arrow} is {arrow.source} which is not a vertex {self.vertex_list}.")
            try:
                target_idx = self.vertex_to_index[arrow.target]
            except:
                raise ValueError(f"Target of {arrow} is {arrow.target} which is not a vertex {self.vertex_list}.")
            self.succ_masks[source_idx] |= (1 << target_idx)
            self.pred_masks[target_idx] |= (1 << source_idx)
            self.adj_masks[source_idx] |= (1 << target_idx)
            self.adj_masks[target_idx] |= (1 << source_idx)

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

    def _sources_of_mask(self, mask: int) -> int:
        """
        Return bitmask of source vertices in the given mask.
        A source vertex is a vertex with no predecessors in the mask.
        """
        sources = 0
        rem_mask = mask
        while rem_mask:
            bit = rem_mask & -rem_mask
            rem_mask &= ~bit
            idx = bit.bit_length() - 1
            if self.pred_masks[idx] & mask == 0:
                sources |= bit
        return sources
    
    @cached_property
    def _sources(self) -> int:
        """
        Return bitmask of vertices that are sources of the entire graph.
        A source vertex is a vertex with no predecessors.
        """
        return self._sources_of_mask(self.vertex_mask)
    
    def sources(self) -> list[Vertex]:
        """
        Return list of vertices that are sources of the graph.
        A source vertex is a vertex with no predecessors in the mask.
        """
        return self._mask_to_vertices(self._sources)
    
    @cached_property
    def _radical_layers(self) -> list[int] | None:
        """
        Returns list of bitmasks representing radical layers.
        Layer 0 consists of all vertices with no predecessors; 
        removing them yields the next layer, and so on.
        """
        remaining_mask = self.vertex_mask
        layers = []
        while remaining_mask:
            sources = self._sources_of_mask(remaining_mask)
            if sources == 0:
                vertices = self._mask_to_vertices(remaining_mask)
                raise ValueError(f"The graph contains a directed cycle in vertices: {[v.label for v in vertices]}")
            layers.append(sources)
            remaining_mask &= ~sources
        return layers
    
    def radical_layers(self) -> list[list[Vertex]]:
        """
        Returns list of vertices representing radical layers.
        Layer 0 consists of all vertices with no predecessors; 
        removing them yields the next layer, and so on.
        """
        layers = []
        if self._radical_layers:
            for layer_mask in self._radical_layers:
                layers.append(self._mask_to_vertices(layer_mask))
        return layers
    
    def compute_radical_subgraphs(self) -> list[list[Vertex]]:
        """ 
        Returns list of connected components of the radical subgraphs as lists of vertices.
        A radical subgraph at i is the full subgraph of all the vertices in radical layers j>=i. 
        """
        rad_layers = self._radical_layers
        if not rad_layers:
            return []
        subgraph_mask = 0
        subgraphs = []
        for layer_mask in reversed(rad_layers):
            new_subgraph_mask = subgraph_mask + layer_mask
            new_conn_comp_masks = self._connected_components_of(new_subgraph_mask)
            new_conn_comp = [self._mask_to_vertices(mask) for mask in new_conn_comp_masks]
            subgraphs += new_conn_comp
            subgraph_mask = new_subgraph_mask
        subgraphs.reverse()
        return subgraphs

    # Socle layers

    def _sinks_of_mask(self, mask: int) -> int:
        """
        Return bitmask of sink vertices in the given mask.
        A sink vertex is one with no successors in the mask.
        """
        sinks = 0
        for bit in self._iterate_over_bits(mask):
            idx = bit.bit_length() - 1
            if self.succ_masks[idx] & mask == 0:
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
        Layer 0 consists of all vertices with no successors; 
        removing them yields the next layer, and so on.
        """
        remaining_mask = self.vertex_mask
        layers = []
        while remaining_mask:
            sinks = self._sinks_of_mask(remaining_mask)
            if sinks == 0:
                vertices = self._mask_to_vertices(remaining_mask)
                raise ValueError(f"The graph contains a directed cycle in vertices: {[v.label for v in vertices]}")
            layers.append(sinks)
            remaining_mask &= ~sinks
        return layers
    
    def socle_layers(self) -> list[list[Vertex]]:
        """
        Returns list of vertices representing socle layers.
        Layer 0 consists of all vertices with no successors; 
        removing them yields the next layer, and so on.
        """
        layers = []
        for layer_mask in self._socle_layers:
            layers.append(self._mask_to_vertices(layer_mask))
        return layers
    
    ### Connected components ###
    
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
            if self.adj_masks[idx] & smaller_mask:
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
                neighbours = mask & self.adj_masks[idx]
                new_neighbours = neighbours & ~component
                component |= new_neighbours
                to_check |= new_neighbours
            components.append(component)
            remaining_mask &= ~component
        return components
    
    @cached_property
    def _connected_components(self) -> list[int]:
        """ Return list of connected components of graph. """
        return self._connected_components_of(self.vertex_mask)
    
    @cached_property
    def connected_components(self) -> list[list[Vertex]]:
        """ Return list of connected components of graph as lists of vertices. """
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
            pred_v = self.pred_masks[v_idx]
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
            succ_v = self.succ_masks[v_idx]
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
            - predecessor-closed iff  
                for every vertex v in S, and for all arrows u -> v in the graph, u is in S.
            - successor-closed iff  
                for every vertex v in S, and for all arrows v -> u in the graph, u is in S.

        Only connected subsets are considered;
            disconnected masks are ignored and left out of the returned lists.

        Returns:
            tuple[list[list[int]],list[list[int]]]:
                - 0 entry = list of predecessor closed connected vertex subsets.
                - 1 entry = list of successor closed connected vertex subsets.
        """        
        full_mask = self.vertex_mask
        pred_masks = self.pred_masks
        succ_masks = self.succ_masks

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
                    if pred_masks[vertex_idx] & (~mask): # check if p_mask outside of mask
                        pred_closed = False
                    if succ_masks[vertex_idx] & (~mask): # check if s_mask outside of mask
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

        A subset S is predecessor-closed iff  
            for every vertex v in S, and for all arrows u -> v in the graph, u is in S.

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

        A subset S is successor-closed iff  
            for every vertex v in S, and for all arrows v -> u, u is in S.

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
    
    import csv

    class BitmaskGraphExamples:
        """
        Stores a collection of example bitmask graphs for testing.

        Each example is identified by a string key and maps to a tuple containing:
            -   A list of vertices (list[Vertex]) or None.
                If None, the vertex set is the set of all sources and targets of arrows
                in the arrow list.
            -   A list of arrows (list[Arrow]) or None.
                If None, the graph consists only of the given vertices with no edges.
        """
        def __init__(self, 
                    examples: dict[str, tuple[list[Vertex] | None, list[Arrow] | None]]):
            self.examples = examples

        def add(self, 
                name: str, 
                vertex_list: list[Vertex] | None = None, 
                arrow_list: list[Arrow] | None = None):
            self.examples[name] = (vertex_list, arrow_list)

        def add_random_graph(self, name: str, num_vertices: int, edge_probability: float):
            """
            Add a random graph with the given number of vertices and edge probability.

            For each pair of vertices (u,v) with u < v, an arrow u -> v is included
            with probability edge_probability independently of all other arrows.
            """
            vertex_list = [Vertex(str(i)) for i in range(num_vertices)]
            label_to_vertex = {v.label : v for v in vertex_list}
            arrow_list = []
            for i in range(num_vertices):
                for j in range(i+1,num_vertices):
                    if random.random() < edge_probability:
                        arrow_list.append(Arrow(label_to_vertex[str(i)], label_to_vertex[str(j)]))
            self.add(name, vertex_list, arrow_list)

        def run(self, 
                print_example_names: bool = True,
                print_times: bool = True,
                print_results: bool = False) -> dict[str, float]:
            """
            Execute all stored example graphs and measure computation times.

            For each example, this method:
                1. Constructs the corresponding BitmaskGraph.
                2. Computes connected subgraphs.
                3. Compute radical layers and radical subgraphs.
                4. Compute socle layers.
                5. Computes connected predecessor and successor closed subgraphs.
            """

            times = {}
            for name, (vertex_list, arrow_list) in self.examples.items():
                if print_example_names:
                    print(f"\n=== BitmaskGraph Example: {name} ===")
                if not vertex_list:
                    if not arrow_list:
                        print("This example has neither a list of vertices, nor a list of arrows.")
                        continue
                    else:
                        sources = set(a.source for a in arrow_list)
                        targets = set(a.target for a in arrow_list)
                        vertex_list = list( sources | targets )

                start_time = time.time()
                bitmask_graph = BitmaskGraph(vertex_list, arrow_list)
                init_time = time.time() - start_time
                times["initialisation"] = init_time

                if print_times:
                    print(f"\nTime to initialise graph: {init_time:.5f} seconds")

                start_time = time.time()
                bitmask_graph._connected_masks
                conn_time = time.time() - start_time
                times["connected_masks"] = conn_time

                if print_times:
                    print(f"Time to compute list of connected masks: {conn_time:.5f} seconds")

                start_time = time.time()
                connected_comps = bitmask_graph.connected_components
                comp_time = time.time() - start_time
                times["connected_components"] = comp_time

                if print_times:
                    print(f"Time to compute list of connected components: {comp_time:.5f} seconds")

                start_time = time.time()
                rad_layers = bitmask_graph.radical_layers()
                rad_time = time.time() - start_time
                times["radical_layers"] = rad_time

                if print_times:
                    print(f"Time to compute radical layers {rad_time:.5f} seconds")

                start_time = time.time()
                rad_subgraphs = bitmask_graph.compute_radical_subgraphs()
                rad_sub_time = time.time() - start_time
                times["radical_subgraphs"] = rad_sub_time

                if print_times:
                    print(f"Time to compute radical subgraphs: {rad_sub_time:.5f} seconds")

                start_time = time.time()
                soc_layers = bitmask_graph.socle_layers()
                soc_time = time.time() - start_time
                times["socle_layers"] = soc_time

                if print_times:
                    print(f"Time to compute socle layers {soc_time:.5f} seconds")

                start_time = time.time()
                closed = bitmask_graph.compute_closed_subsets
                closed_time = time.time() - start_time
                times["connected_closed_masks"] = closed_time

                if print_times:
                    print(f"Time to compute closed subsets: {closed_time:.5f} seconds")

                if print_results:

                    print("\nVertices (with masks):")
                    for vert, mask in bitmask_graph.vertex_to_index.items():
                        print(f" {vert} -> Mask: {mask}")
                    
                    print("\nArrows:")
                    for arrow in bitmask_graph.arrow_list:
                        if arrow.label:
                            print(f" {arrow.label}: {arrow.source.label} -> {arrow.target.label}")
                        else:
                            print(f" {arrow.source.label} -> {arrow.target.label}")

                    print("\nRadical layers")
                    for idx, layer in enumerate(rad_layers):
                        print(f" Radical layer: {idx} = Vertices: {[v.label for v in layer]}")

                    print("\nRadical graphs")
                    for rad_graph in rad_subgraphs:
                        print(f" Vertices: {[v.label for v in rad_graph]}")

                    print("\nSocle layers")
                    for idx, layer in enumerate(soc_layers):
                        print(f" Socle layer: {idx} = Vertices: {[v.label for v in layer]}")

                    print("\nConnected components:")
                    for comp in connected_comps:
                        print(f" Vertices: {[v.label for v in comp]}")
        
                    print("\nSubmodules:")
                    for list_s in closed[0]:
                        for s in list_s:
                            print(f" Vertices: {[v.label for v in s]}")

                    print("\nQuotients:")
                    for list_p in closed[1]:
                        for p in list_p:
                            print(f" Vertices: {[v.label for v in p]}")

            return times
        
        def collect_time_stats(self, 
                               filename: str,
                               num_runs: int,
                               max_num_vertices: int,
                               num_vertices_increment: int,
                               max_edge_probability: float,
                               edge_probability_increment: float):
            example_counter = 0
            with open(filename, "w", newline="") as f:
                writer = csv.writer(f)
                self.add_random_graph("header_temp", 1, 0.0)
                times = self.run(print_example_names=False, print_times=False)
                if not times:
                    raise ValueError("No method timings returned by run(). Check self.run implementation.")
                time_headers = list(self.run(print_times=False).keys())
                headers = ["num_vertices", "edge_probability", "run_index"] + time_headers
                writer.writerow(headers)
                
                for num_vertices in range(num_vertices_increment, max_num_vertices+num_vertices_increment, num_vertices_increment):
                    print(f"\n=== Collecting data for graphs with {num_vertices} vertices. ===")
                    edge_probabilities = []
                    ep = edge_probability_increment
                    while ep < max_edge_probability + edge_probability_increment:
                        edge_probabilities.append(ep)
                        ep += edge_probability_increment

                    for edge_probability in edge_probabilities:
                        print(f" --- Collecting data for edge probability {edge_probability:.4f} --- ")
                        for run_idx in range(num_runs):
                            self.add_random_graph(f"Example {example_counter}",
                                                  num_vertices,
                                                  edge_probability)
                            example_counter += 1
                            time_data = list(self.run(print_example_names=False,print_times=False).values())
                            data = [num_vertices, edge_probability, run_idx] + time_data
                            writer.writerow(data)
                            print(f"    Finished run {run_idx}.    ")
            
    examples = BitmaskGraphExamples({})

    # for i in range(1):
    #     examples.add_random_graph(f"Random graph {i}", 10, 0.4)

    # times = examples.run(verbose=False)

    examples.collect_time_stats("bitmask_graph_time_stats.csv",
                                5,
                                25,
                                5,
                                0.5,
                                0.1)

    