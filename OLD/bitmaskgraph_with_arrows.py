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
    
    def _mask_to_verts(self, mask : int):
        """
        Convert a vertex bitmask to a list of vertex indices.
        """
        if mask & (~self.full_vertex_mask) != 0:
            raise ValueError(f"Invalid mask. {mask} does not correspond to a list of vertices.")
        return [v for v, m in self.vertex_mask.items() if mask & m != 0 ]
    
    @cached_property
    def adj_mask(self):
        adj = {v : {"p" : 0, "s" : 0} for v in self.vertices}
        vertex_mask = self.vertex_mask
        for a in self.arrows:
            # a is an arrow a.source -> a.target
            # so a.source is a predecessor of a.target
            # and a.target is a successor of a.source
            adj[a.source]["s"] |= vertex_mask[a.target]
            adj[a.target]["p"] |= vertex_mask[a.source]
        return adj

    @cached_property
    def pred_mask(self):
        return {mask : adj_mask["p"] for mask, adj_mask in self.adj_mask.items()}

    @cached_property
    def succ_mask(self):
        return {mask : adj_mask["s"] for mask, adj_mask in self.adj_mask.items()}
    
    @cached_property
    def undirected_adj_mask(self):
        return {mask : adj_mask["p"] & adj_mask["s"] for mask, adj_mask in self.adj_mask.items()}
   
    def is_connected(self, mask : int):
        undir_adj = self.undirected_adj_mask

        # mask is connected iff starting at a single vertex can visit every other vertex
        start = mask & (~mask + 1) # first non-zero bit in mask
        visited = start
        remaining = start

        while remaining != 0:
            
            neighbours = 0
            for mask, adj in undir_adj.items():
                if mask & remaining != 0:
                    neighbours &= adj
            
            neighbour &= mask # only keep neighbours inside mask
            new = neighbours & ~visited # new neighbours that have not been visited
            visited |= new # now they have been visited
            remaining = new # 0 iff neighbours are already in mask or have been visited

        return visited == mask
    
    @cached_property
    def compute_closed_subsets(self):
        pred_mask = self.pred_mask
        succ_mask = self.succ_mask
        pred_closed = {}
        succ_closed = {}
        for mask in range(1, self.full_vertex_mask + 1):
            if not self.is_connected(mask):
                continue
            pred_closed[mask] = False
            p_closed = True
            succ_closed[mask] = False
            s_closed = True

            for single_mask in pred_mask.keys():
                if mask & single_mask != 0:
                    p_mask = pred_mask[single_mask]
                    s_mask = succ_mask[single_mask]
                    
                    if p_mask & (~mask): # check if p_mask outside of mask
                        p_closed = False
                    if s_mask & (~mask): # check if s_mask outside of mask
                        s_closed = False
                    
                    if not(p_closed or s_closed):
                        break

            if p_closed:
                pred_closed[mask] = True
            if s_closed:
                succ_closed[mask] = True
        return {"pred_closed" : pred_closed, "succ_closed" : succ_closed}





    # @cached_property
    # def _mask_to_vertex_lists(self):
    #     """
    #     Convert each vertex mask to a list of vertex indices.
    #     """
    #     return {mask : self._mask_to_verts(mask) for mask in range(1, self._mask_all_vertices+1)}
    
    # @cached_property
    # def _mask_vertex_lists(self):
    #     """
    #     Convert each list of vertex indices to a vertex mask.
    #     """
    #     return {v : k for k, v in self._mask_to_vertex_lists.items()}

    
    # ### ARROWS ###
    
    # @cached_property
    # def _mask_arrows(self) -> dict[Arrow,tuple]:
    #     """
    #     Compute bitmask for each arrow and cache.

    #     For V vertices and A arrows: arrows in 0..V+A bits

    #     The bitmask of an arrow has four components:
    #     - Source vertex bit (in 0..V-1 bits)
    #     - Target vertex bit (in 0..V-1 bits)
    #     - Unique arrow bit (in V..A-1 bits)
    #     - Directional bit in position 0 (in V+A bit)
    #         - 0 if source bit > target bit
    #         - 1 if source bit < target bit

    #     Returns:
    #         arrow -> (mask, source mask, target mask, arrow mask)
    #     """
    #     mask_verts = self._mask_vertices
    #     V = self.num_vertices
    #     A = self.num_arrows
    #     VplusA = 1 << V+A
        
    #     mask_arrs = {}

    #     for counter, arrow in enumerate(self.arrows):
    #         if arrow.source == arrow.target:
    #             raise ValueError(f"Invalid arrow {arrow}. Cannot contain loops.")

    #         source_bit = mask_verts[arrow.source]
    #         target_bit = mask_verts[arrow.target]

    #         directional_bit = VplusA if source_bit < target_bit else 0

    #         arrow_bit = 1 << V + counter

    #         full_mask = source_bit | target_bit | arrow_bit | directional_bit
            
    #         mask_arrs[arrow] = (full_mask,source_bit,target_bit,arrow_bit)

    #     return mask_arrs
    
    # @cached_property
    # def _mask_to_arrow_info(self):
    #     """
    #     Compute the source, target and unique arrow bit from a full mask.

    #     full mask -> (source bit, target bit, arrow bit)
    #     """
    #     return {mask : (source,target,arrow) for (mask,source,target,arrow) in self._mask_arrows.values()}
    
    # def _arrow_info_from_mask(self, mask : int):
    #     """
    #     Compute source, target of an arrow from its mask.
    #     """
    #     arrows = []
    #     for arrow_obj, (arrow_mask, _, _, _) in self._mask_arrows.values():
    #         if mask & arrow_mask:
    #             arrows.append(arrow_obj)
    #     return arrows
    
    # ### ADJACENCIES ###

    # @cached_property
    # def _mask_single_vertex_adjacencies(self):
    #     """
    #     Map each single vertex mask to its adjacency info:
    #         mask -> {
    #             "pred_info" : {"verts" : vertex mask, "arrs" : arrow mask}, 
    #             "succ_info" : {"verts" : vertex mask, "arrs" : arrow mask}
    #             }
            
    #     Here:
    #         pred_info: vertices/arrows that point into this vertex
    #         succ_info: vertices/arrows that this vertex points to
    #     """
    #     single_adj = {
    #         mask : {
    #             "pred_info" : {"verts" : 0, "arrs" : 0}, 
    #             "succ_info" : {"verts" : 0, "arrs" : 0}
    #             } 
    #             for mask in self._mask_vertices.values()}
    #     for (full_mask, source_mask, target_mask, _) in self._mask_arrows.values():
    #         # there is an arrow source_mask -> target_mask
    #         # so target_mask is a successor of source_mask
    #         single_adj[source_mask]["succ_info"]["verts"] |= target_mask
    #         single_adj[source_mask]["succ_info"]["arrs"] |= full_mask
    #         # and source_mask is a predecessor of target_mask
    #         single_adj[target_mask]["pred_info"]["verts"] |= source_mask
    #         single_adj[target_mask]["pred_info"]["arrs"] |= full_mask
    #     return single_adj

    # def _update_adj_info(self, list_adjs, list_mask, single_mask, single_adj):
    #     """
    #     Update adjacency info for list_mask in list_adjs using the adjacency of single_mask.

    #     We inherit predecessor/successor info from each single vertex *in* the list
    #     but we DO NOT add any predecessor/successor vertex that is already inside
    #     list_mask.
    #     """
    #     # Skip unless this single vertex is part of the list
    #     if not (list_mask & single_mask):
    #         return

    #     pred_info = single_adj["pred_info"]
    #     succ_info = single_adj["succ_info"]

    #     # Add predecessors not already in the list
    #     if not (list_mask & pred_info["verts"]):
    #         list_adjs[list_mask]["pred_info"]["verts"] |= pred_info["verts"]
    #         list_adjs[list_mask]["pred_info"]["arrs"]  |= pred_info["arrs"]

    #     # Add successors not already in the list
    #     if not (list_mask & succ_info["verts"]):
    #         list_adjs[list_mask]["succ_info"]["verts"] |= succ_info["verts"]
    #         list_adjs[list_mask]["succ_info"]["arrs"]  |= succ_info["arrs"]


    # @cached_property
    # def _mask_vertex_list_adjacencies(self):
    #     """
    #     Map each list of vertices mask to its adjacency info:
    #         mask -> {
    #             "pred_info" : {"verts" : vertex mask, "arrs" : arrow mask}, 
    #             "succ_info" : {"verts" : vertex mask, "arrs" : arrow mask}
    #             }
        
    #     A list of vertices inherits adjacencies from each of its entries.
    #     Note we do not add a pred/succ vertex if it is already inside the list.
    #     """
    #     single_adjs = self._mask_single_vertex_adjacencies
    #     full_mask = self._mask_all_vertices

    #     list_adjs = {
    #         mask : {
    #             "pred_info" : {"verts" : 0, "arrs" : 0}, 
    #             "succ_info" : {"verts" : 0, "arrs" : 0}
    #             } 
    #             for mask in range(1, full_mask + 1)}
        
    #     for list_mask in range(1, full_mask + 1):
    #         for single_mask, single_adj in single_adjs.items():
    #             self._update_adj_info(list_adjs, list_mask, single_mask, single_adj)
    #     return list_adjs
    
    # @cached_property
    # def _DIRECT_mask_vertex_list_adjacencies(self):
    #     """
    #     vertex list mask -> {pred_info : {verts, arrs}, succ_info : (vertices, arrows)}
    #     """
    #     adj = {mask : {"pred_info" : {"verts" : 0, "arrs" : 0}, "succ_info" : {"verts" : 0, "arrs" : 0}} for mask in range(1,self._mask_all_vertices+1)}
    #     for verts in range(1,self._mask_all_vertices+1):
    #         for (full_mask, source_mask, target_mask, _) in self._mask_arrows.values():
    #             if verts & source_mask != 0 and verts & target_mask == 0:
    #                 adj[verts]["succ_info"]["verts"] |= target_mask
    #                 adj[verts]["succ_info"]["arrs"] |= full_mask
    #             if verts & target_mask != 0 and verts & source_mask == 0:
    #                 adj[verts]["pred_info"]["verts"] |= source_mask
    #                 adj[verts]["pred_info"]["arrs"] |= full_mask
    #     return adj

    # ### CONNECTIVITY ###

    # @cached_property
    # def _connected_vertex_lists(self):
    #     """
        
    #     """
    #     single_adjs = self._mask_single_vertex_adjacencies
    #     list_adjs = self._mask_vertex_list_adjacencies
    #     connected = {}
    #     for n in range(1, self._mask_all_vertices + 1):
    #         connected[n] = False
    #         start = n & (~n + 1)
    #         visited = start
    #         remaining = start

    #         while remaining != 0:
                
    #             neighbours = 0
    #             for single_mask, single_adj in single_adjs.items():
    #                 if remaining & single_mask != 0:
    #                     neighbours |= single_adj["pred_info"]["verts"]
    #                     neighbours |= single_adj["succ_info"]["verts"]
    #             neighbours &= n # only keep neighbours inside n

    #             new = neighbours & ~visited # new neighbours that have not been visited
    #             visited |= new
    #             remaining = new
    #         connected[n] = (visited == n)
    #     return connected

    # ### CONNECTED SUBGRAPHS/QUOTIENTS ###
        
    # @cached_property
    # def _mask_vertex_list_closed_adjacencies(self):
    #     mask_adj = self._mask_vertex_list_adjacencies
    #     conn_lists = self._connected_vertex_lists
    #     pred_closed = []
    #     succ_closed = []
    #     for mask, info in mask_adj.items():
    #         pred_verts = info["pred_info"]["verts"]
    #         succ_verts = info["succ_info"]["verts"]
    #         if pred_verts == 0 and conn_lists[mask]:
    #             pred_closed.append(mask)
    #         if succ_verts == 0 and conn_lists[mask]:
    #             succ_closed.append(mask)
    #     return {"pred_closed" : pred_closed, "succ_closed" : succ_closed}




    # ### TODO: which of these functions are useful? ###

                




    # @cached_property
    # def _mask_vertex_list_predecessors(self):
    #     """
    #     vertex list mask -> {verts : predecessor vertices, arrs : predecessor arrows}
    #     """
    #     return {k : v["pred_info"] for k, v in self._mask_vertex_list_adjacencies.items()} 
    
    # @cached_property
    # def _mask_vertex_list_successors(self):
    #     """
    #     vertex list mask -> {verts : successor vertices, arrs : successor arrows}
    #     """
    #     return {k : v["succ_info"] for k, v in self._mask_vertex_list_adjacencies.items()} 


    # @cached_property
    # def _mask_source_vertices(self):
    #     """
    #     Return mask of source vertices.
    #     """
    #     source = 0
    #     mask_pred = self._mask_vertex_list_predecessors
    #     for info in mask_pred.values():
    #         pred_verts = info["verts"]
    #         if pred_verts == 0:
    #             source |= pred_verts
    #     return source
    
    # @cached_property
    # def _mask_sink_vertices(self):
    #     """
    #     Return mask of sink vertices.
    #     """
    #     sink = 0
    #     mask_succ = self._mask_vertex_list_successors
    #     for info in mask_succ.values():
    #         succ_verts = info["verts"]
    #         if succ_verts == 0:
    #             source |= succ_verts
    #     return source


    # @cached_property
    # def _mask_predecessors(self):
    #     """
    #     vertex mask -> {predecessor vertices, predecessor arrows}
    #     """
    #     return {k : v["pred_info"] for k, v in self._mask_vertex_adjacencies.items()} 
    
    # @cached_property
    # def _mask_successors(self):
    #     """
    #     vertex mask -> {successor vertices, successor arrows}
    #     """
    #     return {k : v["succ_info"] for k, v in self._mask_vertex_adjacencies.items()} 
    
    # @cached_property
    # def _mask_source_vertices(self):
    #     """
    #     Compute bitmask corresponding to source vertices.
    #     """
    #     sources = 0
    #     for mask, pred_info in self._mask_predecessors.items():
    #         if pred_info["verts"] == 0:
    #             sources |= mask
    #     return sources
    
    # @cached_property
    # def compute_radical_layers(self):
    #     """
    #     Use bitmask to compute radical layers.
    #     """
    #     mask_succ = self._mask_successors
    #     sources = self._mask_source_vertices
    #     rad_layer = {0 : sources}
    #     all_verts = self._mask_all_vertices
    #     test_verts = all_verts & (~sources)
    #     previous_layer = 0
    #     while test_verts != 0:
    #         previous_verts = rad_layer[previous_layer]

    #         return


    # # @cached_property
    # # def _mask_sink_vertices(self):
    # #     """
    # #     Compute bitmask corresponding to sink vertices.
    # #     """
    # #     sinks = 0
    # #     for mask, succ_info in self._mask_successors.items():
    # #         if succ_info["verts"] == 0:
    # #             sinks |= mask
    # #     return sinks
        
    # def _mask_adjacency_of_vertex(self, mask_vert : int):
    #     """
    #     Compute adjacency for a specific vertex mask.
    #     """
    #     if mask_vert & self._mask_all_vertices == 0:
    #         raise ValueError(f"Invalid mask_vert.")
    #     return self._mask_vertex_list_adjacencies[mask_vert]
    
    # def _mask_predecessors_of_vertex(self, mask_vert : int):
    #     """
    #     Compute predecessor information for a specific vertex mask.
    #     """
    #     return self._mask_vertex_list_predecessors[mask_vert]
        
    # def _mask_successesors_of_vertex(self, mask_vert : int):
    #     """
    #     Computer successor information for a specific vertex mask.
    #     """
    #     return self._mask_vertex_list_successors[mask_vert]
    


    # # TODO: can this be done in one go?
    # def _mask_anc_info(self, mask_vert : int):
    #     """
    #     Compute bitmask corresponding to all ancestor vertices for a specific vertex mask.
    #     """
    #     mask_adj = self._mask_vertex_adjacencies
    #     mask_anc_verts = mask_vert
    #     mask_anc_arrs = 0
    #     mask_test = mask_vert
    #     mask_verts = self._mask_vertices
    #     while mask_test != 0:
    #         for mask in mask_verts.values():
    #             if mask_test & mask != 0:
    #                 mask_pred = mask_adj[mask]["pred_info"]
    #                 mask_anc_verts |= mask_pred["verts"]
    #                 mask_anc_arrs |= mask_pred["arrs"]
    #                 mask_test = mask_test & (~ mask) | mask_pred["verts"]
    #     return {"verts" : mask_anc_verts, "arrs" : mask_anc_arrs}
    
    # def _mask_des_info(self, mask_vert : int):
    #     """
    #     Compute bitmask corresponding to all descendant vertices for a specific vertex mask.
    #     """
    #     mask_adjacency = self._mask_vertex_adjacencies
    #     mask_des_verts = mask_vert
    #     mask_des_arrs = 0
    #     mask_test = mask_vert
    #     mask_verts = self._mask_vertices
    #     while mask_test != 0:
    #         for mask in mask_verts.values():
    #             if mask_test & mask != 0:
    #                 mask_succ = mask_adjacency[mask]["succ_info"]
    #                 mask_des_verts |= mask_succ["verts"]
    #                 mask_des_arrs |= mask_succ["arrs"]
    #                 mask_test = mask_test & (~ mask) | mask_succ["verts"]
    #     return {"verts" : mask_des_verts, "arrs" : mask_des_arrs}
    
    # def _mask_anc_verts(self, mask_vert : int):
    #     return self._mask_anc_info(mask_vert)["verts"]
    
    # def _mask_anc_arrs(self, mask_vert : int):
    #     return self._mask_anc_info(mask_vert)["arrs"]
    
    # def _mask_des_verts(self, mask_vert : int):
    #     return self._mask_des_info(mask_vert)["verts"]
    
    # def _mask_des_arrs(self, mask_vert : int):
    #     return self._mask_des_info(mask_vert)["arrs"]
    
    # def _mask_anc_verts_list(self, mask_verts : int):
    #     return sum(self._mask_anc_verts(v) for v in self._mask_vertices if v & mask_verts != 0)
    
    # def _mask_des_verts_list(self, mask_verts : int):
    #     return sum(self._mask_des_verts(v) for v in self._mask_vertices if v & mask_verts != 0)

    # def _mask_anc_verts_comp(self, mask_vert_list : int):
    #     comp = []
    #     mask_verts = self._mask_vertices
    #     mask_source = self._mask_source_vertices
    #     mask_rem_verts = mask_vert_list
    #     mask_result = 0
    #     while mask_rem_verts != 0:
    #         for mask in mask_verts.values():
    #             if mask & mask_rem_verts != 0:
    #                 mask_result |= self._mask_anc_verts(mask)
    #                 mask_rem_verts &= (~mask)
    #                 mask_s = mask_source & mask_result
    #                 if self._mask_des_verts(mask_s) & mask_rem_verts == 0:
    #                     comp.append(mask_result)
    #                     mask_result = 0
    #     return comp
    
    # def _mask_des_verts_comp(self, mask_vert_list : int):
    #     comp = []
    #     mask_verts = self._mask_vertices
    #     mask_sink = self._mask_sink_vertices
    #     mask_rem_verts = mask_vert_list
    #     mask_result = 0
    #     while mask_rem_verts != 0:
    #         for mask in mask_verts.values():
    #             if mask & mask_rem_verts != 0:
    #                 mask_result |= self._mask_des_verts(mask)
    #                 mask_rem_verts &= (~mask)
    #                 mask_s = mask_sink & mask_result
    #                 if self._mask_anc_verts(mask_s) & mask_rem_verts == 0:
    #                     comp.append(mask_result)
    #                     mask_result = 0
    #     return comp
    
    # @cached_property
    # def is_connected(self):
    #     return len(self._mask_des_verts_comp(self._mask_all_vertices)) == 1

    # @cached_property
    # def _mask_anc_closed(self):
    #     """
    #     Find all ancestor closed subsets of vertices.
    #     """
    #     result = []
    #     mask_verts = self._mask_vertices
    #     all_verts = self._mask_all_vertices
    #     for test_mod in range(1,all_verts+1):
    #         left_to_check = test_mod
    #         for mask_vert in mask_verts.values():
    #             if mask_vert & test_mod == mask_vert:
    #                 mask_anc_vert = self._mask_anc_verts(mask_vert)
    #                 if mask_anc_vert & test_mod == mask_anc_vert:
    #                     left_to_check &= (~mask_anc_vert)
    #                     if left_to_check == 0:
    #                         comp = self._mask_anc_verts_comp(test_mod)
    #                         if len(comp) == 1:
    #                             result.append(test_mod)
    #                         break
    #                 else:
    #                     break           
    #     return result

    # @cached_property
    # def _mask_des_closed(self):
    #     result = []
    #     mask_verts = self._mask_vertices
    #     all_verts = self._mask_all_vertices
    #     for test_mod in range(1,all_verts+1):            
    #         left_to_check = test_mod
    #         for mask_vert in mask_verts.values():
    #             if mask_vert & test_mod == mask_vert:
    #                 mask_des_vert = self._mask_des_verts(mask_vert)
    #                 if mask_des_vert & test_mod == mask_des_vert:
    #                     left_to_check &= (~mask_des_vert)
    #                     if left_to_check == 0:
    #                         comp = self._mask_des_verts_comp(test_mod)
    #                         if len(comp) == 1:
    #                             result.append(test_mod)
    #                         break
    #                 else:
    #                     break
    #     return result
    
    # @cached_property
    # def anc_closed(self):
    #     results = self._mask_anc_closed
    #     anc_closed = [[]]
    #     for result in results:
    #         anc_closed.append(self._mask_to_verts(result))
    #     return anc_closed
    
    # @cached_property
    # def des_closed(self):
    #     results = self._mask_des_closed
    #     des_closed = [[]]
    #     for result in results:
    #         des_closed.append(self._mask_to_verts(result))
    #     return des_closed
    
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

            mask_to_lists = graph._mask_to_vertex_lists
            print("\nVertex lists (with masks):")
            for mask, list in mask_to_lists.items():
                print(f"  {list} -> mask: {mask}")

            print("\nArrows:")
            for arr, mask in graph._mask_arrows.items():
                print(f"  {arr} -> mask: {mask}")

            print("\n Connected subgraphs:")
            for n, is_conn in graph._connected_vertex_lists.items():
                if is_conn:
                    print(f" {n} : {mask_to_lists[n]}")

            adj_closed = graph._mask_vertex_list_closed_adjacencies
            print("\nSubmodules:")
            for s in adj_closed["succ_closed"]:
                print(f"  {s} : {mask_to_lists[s]}")

            print("\nQuotients:")
            for q in adj_closed["pred_closed"]:
                print(f"  {q} : {mask_to_lists[q]}")

        

if __name__ == "__main__":

    examples = Examples({})

    examples.add("Type A", [Arrow(0,1), Arrow(1,2), Arrow(2,3)])

    examples.add("Diamond", [Arrow(0,1), Arrow(0,2), Arrow(1,3), Arrow(2,3)])

    examples.run()

    a = [Arrow(0,1), Arrow(1,2), Arrow(2,3)]
    b = BitmaskGraph(a)
    print(b.vertices)



