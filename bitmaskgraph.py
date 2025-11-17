from typing import Optional
from functools import cached_property

class Arrow:
    def __init__(self, source : int, target : int, label : Optional[str]=None):
        self.source = source
        self.target = target
        self.label = label

    def __eq__(self, other : Arrow):
        return (self.source == other.source) and (self.target == other.target) and (self.label == other.label)

    def __repr__(self):
        return f"{self.label}: {self.source}->{self.target}"
    
    def __hash__(self):
        return hash((self.source, self.target, self.label))
    
class Bitmask:
    """
    Contains operations to apply to a generic bitmask.
    """
    def __init__(self,
                 info : list):
        self.info = info

    def __len__(self):
        """
        Length of the bitmask is the number of pieces of information to encode.
        """
        return len(self.info)

    @cached_property
    def mask_info(self):
        """
        Mask the information in len(info) number of bits.
        """
        return {x : 1 << idx for idx, x in enumerate(self.info)}
    
    @cached_property
    def mask_all_info(self) -> int:
        """
        Compute bitmask of full list of info.
        """
        return (1 << len(self)) - 1
    
    def mask_to_info(self, mask : int):
        """
        Convert each mask to the corresponding sublist of info.
        """
        if mask & (~ self.mask_all_info) != 0:
            raise ValueError(f"Invalid mask. {mask} does not correspond to a sublist of {self.info}.")
        return [x for x, m in self.mask_info.items() if mask & m != 0 ]
    
    @cached_property
    def mask_to_info_sublists(self):
        """
        Convert each info mask to a sublist of info.
        """
        return {mask : self.mask_to_info(mask) for mask in range(1, self.mask_all_info+1)}
    
    @cached_property
    def mask_info_sublists(self):
        """
        Convert each sublist of info into an info mask.
        """
        return {v : k for k, v in self.mask_to_info_sublists.items()}



class BitmaskGraph:
    """
    Represents a graph using bitmasks.

    For V vertices and A arrows.

    Layout:
    - Bits 0..V-1 : vertex bits
    - Bits V..A-1 : arrow bits
    - Bit V+A : arrow direction bit (1 if source<target, 0 if target<source)

    Attributes:

    """


    def __init__(self,
                 arrows : list[Arrow]):
        self.arrows = arrows
        self.num_arrows = len(arrows)
        self.vertices = list(set(a.source for a in self.arrows) | set(a.target for a in self.arrows))
        self.vertex_to_index = {v : idx for idx, v in enumerate(self.vertices)}
        self.num_vertices = len(self.vertices)

    ### VERTICES ###

    @cached_property
    def vertex_mask(self) -> dict[int,int]:
        """
        Create bitmasks for each vertex and cache.
        For V vertices: vertex mask in bits 0..V-1.
        """
        return {v : 1 << idx for v, idx in self.vertex_to_index.items()}

    @cached_property
    def full_vertex_mask(self) -> int:
        """
        Compute bitmask of all vertices.
        """
        return (1 << self.num_vertices) - 1
    
    def mask_to_verts(self, mask : int):
        """
        Convert a vertex bitmask to a list of vertex indices.
        """
        if mask & (~self.full_vertex_mask) != 0:
            raise ValueError(f"Invalid mask. {mask} does not correspond to a list of vertices.")
        return [v for v, m in self.vertex_mask.items() if mask & m != 0]
    
    @cached_property
    def adj_mask(self):
        vertex_mask = self.vertex_mask
        adj = {vertex_mask[v] : {"p" : 0, "s" : 0} for v in self.vertices}
        vertex_mask = self.vertex_mask
        for a in self.arrows:
            # a is an arrow a.source -> a.target
            # so a.source is a predecessor of a.target
            # and a.target is a successor of a.source
            adj[vertex_mask[a.source]]["s"] |= vertex_mask[a.target]
            adj[vertex_mask[a.target]]["p"] |= vertex_mask[a.source]
        return adj

    @cached_property
    def pred_mask(self):
        return {mask : adj_mask["p"] for mask, adj_mask in self.adj_mask.items()}

    @cached_property
    def succ_mask(self):
        return {mask : adj_mask["s"] for mask, adj_mask in self.adj_mask.items()}
    
    @cached_property
    def undirected_adj_mask(self):
        return {mask : adj_mask["p"] | adj_mask["s"] for mask, adj_mask in self.adj_mask.items()}
   
    def is_connected(self, mask : int):
        undir_adj = self.undirected_adj_mask

        # mask is connected iff starting at a single vertex can visit every other vertex
        start = mask & (~mask + 1) # first non-zero bit in mask
        visited = start
        remaining = start

        while remaining != 0:
            
            neighbours = 0
            for single_mask, adj in undir_adj.items():
                if single_mask & remaining != 0:
                    neighbours |= adj
            
            neighbours &= mask # only keep neighbours inside mask
            new = neighbours & ~visited # new neighbours that have not been visited
            visited |= new # now they have been visited
            remaining = new # 0 iff neighbours are already in mask or have been visited

        return visited == mask
    
    @cached_property
    def compute_closed_subsets(self):

        vert_mask = self.vertex_mask
        full_mask = self.full_vertex_mask
        pred_mask = self.pred_mask
        succ_mask = self.succ_mask
        pred_closed = {}
        succ_closed = {}

        for mask in range(1, full_mask + 1):
            if not self.is_connected(mask):
                continue

            p_closed = True
            s_closed = True

            for single_mask in vert_mask.values():
                if mask & single_mask != 0:
                    p_mask = pred_mask[single_mask]
                    s_mask = succ_mask[single_mask]
                    
                    if p_closed and p_mask & (~mask): # check if p_mask outside of mask
                        p_closed = False
                    if s_closed and s_mask & (~mask): # check if s_mask outside of mask
                        s_closed = False
                    
                    if not(p_closed or s_closed):
                        break

            pred_closed[mask] = p_closed
            succ_closed[mask] = s_closed

        return {"pred_closed" : pred_closed, "succ_closed" : succ_closed}

class Examples:
    """
    Class to store examples of bitmask graphs.
    """
    def __init__(self, 
                 examples: dict[str, list[Arrow]]):
        self.examples = examples

    def add(self, name: str, arrows: list[Arrow]):
        self.examples[name] = arrows

    def run(self):
        for name, arrows in self.examples.items():
            print(f"\n=== BitmaskGraph Example: {name} ===")
            graph = BitmaskGraph(arrows)
            
            print("\nVertices (with masks):")
            for vert, mask in graph.vertex_mask.items():
                print(f"  {vert} -> mask: {mask}")

            print("\nVertex masks (with lists):")
            for list in range(1, graph.full_vertex_mask + 1):
                print(f"  {list} -> vertices: {graph.mask_to_verts(list)}")

            closed_subsets = graph.compute_closed_subsets

            print("\nSubmodules:")
            for s, is_closed in closed_subsets["succ_closed"].items():
                if is_closed:
                    print(f"  {s} : {graph.mask_to_verts(s)}")

            print("\nQuotients:")
            for q, is_closed in closed_subsets["pred_closed"].items():
                if is_closed:
                    print(f"  {q} : {graph.mask_to_verts(q)}")

        

if __name__ == "__main__":

    examples = Examples({})

    examples.add("Type A", [Arrow(0,1), Arrow(1,2), Arrow(2,3)])

    examples.add("Diamond", [Arrow(0,1), Arrow(0,2), Arrow(1,3), Arrow(2,3)])

    examples.run()



