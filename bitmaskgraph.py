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

    @cached_property
    def _mask_vertices(self) -> dict[int,int]:
        """
        Create bitmasks for each vertex and cache.
        For V vertices: vertex mask in bits 0..V-1.
        """
        return {v : 1 << idx for v, idx in self.vertex_to_index.items()}

    @cached_property
    def _mask_all_vertices(self) -> int:
        """
        Compute bitmask of all vertices.
        """
        return sum(self._mask_vertices.values())

    def _mask_to_verts(self, mask : int):
        """
        Convert vertex bitmask to a list of the vertex indices.
        """
        if mask & (~self._mask_all_vertices) != 0:
            raise ValueError(f"Invalid mask. {mask} does not correspond to a list of vertices.")
        return [v for v, m in self._mask_vertices.items() if mask & m != 0 ]
    
    @cached_property
    def _mask_arrows(self) -> dict[Arrow,tuple]:
        """
        Compute bitmask for each arrow and cache.

        For V vertices and A arrows: arrows in 0..V+A bits

        The bitmask of an arrow has four components:
        - Source vertex bit (in 0..V-1 bits)
        - Target vertex bit (in 0..V-1 bits)
        - Unique arrow bit (in V..A-1 bits)
        - Directional bit in position 0 (in V+A bit)
            - 0 if source bit > target bit
            - 1 if source bit < target bit

        Returns:
            arrow -> (mask, source mask, target mask, arrow mask)
        """
        mask_verts = self._mask_vertices
        V = self.num_vertices
        A = self.num_arrows
        VplusA = 1 << V+A
        
        mask_arrs = {}

        for counter, arrow in enumerate(self.arrows):
            if arrow.source == arrow.target:
                raise ValueError(f"Invalid arrow {arrow}. Cannot contain loops.")

            source_bit = mask_verts[arrow.source]
            target_bit = mask_verts[arrow.target]

            directional_bit = VplusA if source_bit < target_bit else 0

            arrow_bit = 1 << V + counter

            full_mask = source_bit | target_bit | arrow_bit | directional_bit
            
            mask_arrs[arrow] = (full_mask,source_bit,target_bit,arrow_bit)

        return mask_arrs
    
    @cached_property
    def _mask_to_arrow_info(self):
        """
        Compute the source, target and unique arrow bit from a full mask.

        full mask -> (source, target, arrow)
        """
        return {mask : (source,target,arrow) for (mask,source,target,arrow) in self._mask_arrows.values()}
    
    def _arrow_info_from_mask(self, mask : int):
        """
        Compute source, target of an arrow from its mask.
        """
        arrows = []
        for arrow_obj, (arrow_mask, _, _, _) in self._mask_arrows.values():
            if mask & arrow_mask:
                arrows.append(arrow_obj)
        return arrows
    
    @cached_property
    def _mask_vertex_adj(self):
        """
        vertex mask -> {pred_info : {verts, arrs}, succ_info : (vertices, arrows)}
        """
        adj = {mask : {"pred_info" : {"verts" : 0, "arrs" : 0}, "succ_info" : {"verts" : 0, "arrs" : 0}} for mask in self._mask_vertices.values()}
        for (full_mask, source_mask, target_mask, _) in self._mask_arrows.values():
            # there is an arrow source_mask -> target_mask
            # so target_mask is a successor of source_mask
            # and source_mask is a predecessor of target_mask
            
            adj[source_mask]["succ_info"]["verts"] |= target_mask
            adj[source_mask]["succ_info"]["arrs"] |= full_mask

            adj[target_mask]["pred_info"]["verts"] |= source_mask
            adj[target_mask]["pred_info"]["arrs"] |= full_mask
        
        return adj
    
    @cached_property
    def _mask_source_vertices(self):
        """
        Compute bitmask corresponding to source vertices.
        """
        sources = 0
        for mask, adj_info in self._mask_vertex_adj.items():
            if adj_info["pred_info"]["verts"] == 0:
                sources |= mask
        return sources
    
    @cached_property
    def _mask_sink_vertices(self):
        """
        Compute bitmask corresponding to sink vertices.
        """
        sinks = 0
        for mask, adj_info in self._mask_vertex_adj.items():
            if adj_info["succ_info"]["verts"] == 0:
                sinks |= mask
        return sinks
        
    def _mask_adj_info(self, mask_vert : int):
        """
        Compute adjacency for a specific vertex mask.
        """
        if mask_vert & self._mask_all_vertices == 0:
            raise ValueError(f"Invalid mask_vert.")
        return self._mask_vertex_adj[mask_vert]
    
    def _mask_pred_info(self, mask_vert : int):
        """
        Compute predecessor information for a specific vertex mask.
        """
        return self._mask_adj_info(mask_vert)["pred_info"]
    
    def _mask_pred_verts(self, mask_vert : int):
        """
        Compute predecessor vertices for a specific vertex mask.
        """
        return self._mask_pred_info(mask_vert)["verts"]
    
    def _mask_pred_arrs(self, mask_vert : int):
        return self._mask_pred_info(mask_vert)["arrs"]
    
    def _mask_succ_info(self, mask_vert : int):
        return self._mask_adj_info(mask_vert)["succ_info"]
    
    def _mask_succ_verts(self, mask_vert : int):
        return self._mask_succ_info(mask_vert)["verts"]
    
    def _mask_succ_arrs(self, mask_vert : int):
        return self._mask_succ_info(mask_vert)["arrs"]

    # TODO: can this be done in one go?
    def _mask_anc_info(self, mask_vert : int):
        """
        Compute bitmask corresponding to all ancestor vertices for a specific vertex mask.
        """
        mask_adj = self._mask_vertex_adj
        mask_anc_verts = mask_vert
        mask_anc_arrs = 0
        mask_test = mask_vert
        mask_verts = self._mask_vertices
        while mask_test != 0:
            for mask in mask_verts.values():
                if mask_test & mask != 0:
                    mask_pred = mask_adj[mask]["pred_info"]
                    mask_anc_verts |= mask_pred["verts"]
                    mask_anc_arrs |= mask_pred["arrs"]
                    mask_test = mask_test & (~ mask) | mask_pred["verts"]
        return {"verts" : mask_anc_verts, "arrs" : mask_anc_arrs}
    
    def _mask_des_info(self, mask_vert : int):
        """
        Compute bitmask corresponding to all descendant vertices for a specific vertex mask.
        """
        mask_adj = self._mask_vertex_adj
        mask_des_verts = mask_vert
        mask_des_arrs = 0
        mask_test = mask_vert
        mask_verts = self._mask_vertices
        while mask_test != 0:
            for mask in mask_verts.values():
                if mask_test & mask != 0:
                    mask_succ = mask_adj[mask]["succ_info"]
                    mask_des_verts |= mask_succ["verts"]
                    mask_des_arrs |= mask_succ["arrs"]
                    mask_test = mask_test & (~ mask) | mask_succ["verts"]
        return {"verts" : mask_des_verts, "arrs" : mask_des_arrs}
    
    def _mask_anc_verts(self, mask_vert : int):
        return self._mask_anc_info(mask_vert)["verts"]
    
    def _mask_anc_arrs(self, mask_vert : int):
        return self._mask_anc_info(mask_vert)["arrs"]
    
    def _mask_des_verts(self, mask_vert : int):
        return self._mask_des_info(mask_vert)["verts"]
    
    def _mask_des_arrs(self, mask_vert : int):
        return self._mask_des_info(mask_vert)["arrs"]
    
    def _mask_anc_verts_list(self, mask_verts : int):
        return sum(self._mask_anc_verts(v) for v in self._mask_vertices if v & mask_verts != 0)
    
    def _mask_des_verts_list(self, mask_verts : int):
        return sum(self._mask_des_verts(v) for v in self._mask_vertices if v & mask_verts != 0)

    def _mask_anc_verts_comp(self, mask_vert_list : int):
        comp = []
        mask_verts = self._mask_vertices
        mask_source = self._mask_source_vertices
        mask_rem_verts = mask_vert_list
        mask_result = 0
        while mask_rem_verts != 0:
            for mask in mask_verts.values():
                if mask & mask_rem_verts != 0:
                    mask_result |= self._mask_anc_verts(mask)
                    mask_rem_verts &= (~mask)
                    mask_s = mask_source & mask_result
                    if self._mask_des_verts(mask_s) & mask_rem_verts == 0:
                        comp.append(mask_result)
                        mask_result = 0
        return comp
    
    def _mask_des_verts_comp(self, mask_vert_list : int):
        comp = []
        mask_verts = self._mask_vertices
        mask_sink = self._mask_sink_vertices
        mask_rem_verts = mask_vert_list
        mask_result = 0
        while mask_rem_verts != 0:
            for mask in mask_verts.values():
                if mask & mask_rem_verts != 0:
                    mask_result |= self._mask_des_verts(mask)
                    mask_rem_verts &= (~mask)
                    mask_s = mask_sink & mask_result
                    if self._mask_anc_verts(mask_s) & mask_rem_verts == 0:
                        comp.append(mask_result)
                        mask_result = 0
        return comp
    
    @cached_property
    def is_connected(self):
        return len(self._mask_des_verts_comp(self._mask_all_vertices)) == 1

    @cached_property
    def _mask_anc_closed(self):
        """
        Find all ancestor closed subsets of vertices.
        """
        result = []
        mask_verts = self._mask_vertices
        all_verts = self._mask_all_vertices
        for test_mod in range(1,all_verts+1):
            left_to_check = test_mod
            for mask_vert in mask_verts.values():
                if mask_vert & test_mod == mask_vert:
                    mask_anc_vert = self._mask_anc_verts(mask_vert)
                    if mask_anc_vert & test_mod == mask_anc_vert:
                        left_to_check &= (~mask_anc_vert)
                        if left_to_check == 0:
                            comp = self._mask_anc_verts_comp(test_mod)
                            if len(comp) == 1:
                                result.append(test_mod)
                            break
                    else:
                        break           
        return result

    @cached_property
    def _mask_des_closed(self):
        result = []
        mask_verts = self._mask_vertices
        all_verts = self._mask_all_vertices
        for test_mod in range(1,all_verts+1):            
            left_to_check = test_mod
            for mask_vert in mask_verts.values():
                if mask_vert & test_mod == mask_vert:
                    mask_des_vert = self._mask_des_verts(mask_vert)
                    if mask_des_vert & test_mod == mask_des_vert:
                        left_to_check &= (~mask_des_vert)
                        if left_to_check == 0:
                            comp = self._mask_des_verts_comp(test_mod)
                            if len(comp) == 1:
                                result.append(test_mod)
                            break
                    else:
                        break
        return result
    
    @cached_property
    def anc_closed(self):
        results = self._mask_anc_closed
        anc_closed = [[]]
        for result in results:
            anc_closed.append(self._mask_to_verts(result))
        return anc_closed
    
    @cached_property
    def des_closed(self):
        results = self._mask_des_closed
        des_closed = [[]]
        for result in results:
            des_closed.append(self._mask_to_verts(result))
        return des_closed
    
    # TODO: this is not working
    # @cached_property
    # def anc_closed_from_des(self):
    #     des_results = self._mask_des_closed
    #     mask_all_verts = self._mask_all_verts
    #     results = [self._mask_to_verts(mask_all_verts)]
    #     for result in des_results:
    #         new_comp = self._mask_des_verts_comp(result)
    #         comp = []
    #         for summand in new_comp:
    #             new_summand = mask_all_verts & (~summand)
    #             comp += self._mask_to_verts(new_summand)
    #         results.append(comp)
    #     return results


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
            for vert, mask in graph._mask_vertices.items():
                print(f"  {vert} -> mask: {mask}")
                print("adj info", graph._mask_adj_info(mask))
                print("anc info", graph._mask_anc_verts(mask))

            print("\nArrows:")
            for arr, mask in graph._mask_arrows.items():
                print(f"  {arr} -> mask: {mask}")



            des_closed = graph.des_closed
            print("\nDescent-closed submodules:")
            for s in des_closed:
                print(f"  {s}")

            anc_closed = graph.anc_closed
            print("\nAncester-closed quotients:")
            for q in anc_closed:
                print(f"  {q}")

        

if __name__ == "__main__":

    examples = Examples({})

    examples.add("Type A", [Arrow(0,1), Arrow(1,2), Arrow(2,3)])

    examples.add("Diamond", [Arrow(0,1), Arrow(0,2), Arrow(1,3), Arrow(2,3)])

    examples.run()

    a = [Arrow(0,1), Arrow(1,2), Arrow(2,3)]
    b = BitmaskGraph(a)
    print(b.vertices)



