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
        return f"Arrow({self.source}->{self.target})"
    
    def __hash__(self):
        return hash((self.source, self.target, self.label))
    

class BitmaskGraph:
    """
    Represents arrows and vertices using bitmasks.

    Attributes:

    """


    def __init__(self,
                 arrows : list[Arrow]):
        self.arrows = arrows
        self.vertex_memory = len(arrows) + 1

    @cached_property
    def _mask_info(self):
        mask_verts = {}
        mask_arrs = {}
        vert_counter = len(self.arrows) + 1
        arr_counter = 1
        for arrow in self.arrows:
            if arrow.source == arrow.target:
                raise ValueError(f"Invalid arrow {arrow}. Cannot contain loops.")
            
            if arrow.source not in mask_verts:
                mask_verts[arrow.source] = 1 << vert_counter
                vert_counter += 1

            if arrow.target not in mask_verts:
                mask_verts[arrow.target] = 1 << vert_counter
                vert_counter += 1

            if arrow not in mask_arrs:
                direction = 0
                if mask_verts[arrow.source] < mask_verts[arrow.target]:
                    direction = 1 << 0
                new_arr_mask = direction | mask_verts[arrow.source] | mask_verts[arrow.target] | (1 << arr_counter)
                mask_arrs[arrow] = new_arr_mask
                arr_counter += 1
        return {"mask_verts" : mask_verts, "mask_arrs" : mask_arrs}
    
    @cached_property
    def _mask_verts(self):
        return self._mask_info["mask_verts"]
    
    def _mask_to_verts(self, mask : int):
        return [v for v, m in self._mask_verts.items() if mask & m != 0 ]
    
    @cached_property
    def _mask_all_verts(self):
        return sum(2**n for n in self._vertices_mask_id)
    
    @cached_property
    def _vertices_mask_id(self):
        return list(range(len(self.arrows)+1, len(self.arrows)+len(self._mask_verts)+1))

    @cached_property
    def _mask_arrs(self):
        return self._mask_info["mask_arrs"]
    
    def _arr_info_from_mask(self, mask_arr : int):
        mask_arrs = self._mask_arrs
        if mask_arr not in mask_arrs.values():
            raise ValueError(f"Invalid mask_arr. {mask_arr} should be in {mask_arrs.values()}.")
        mask_verts = self._mask_verts
        mask_source_target = []
        for mask_vert in mask_verts.values():
            if mask_arr & mask_vert != 0:
                mask_source_target.append(mask_vert)
        mask_min_vert = min(mask for mask in mask_source_target)
        mask_max_vert = max(mask for mask in mask_source_target)
        
        mask_dir = mask_arr & 1
        if mask_dir == 0:
            # direction large to small
            source = mask_max_vert
            target = mask_min_vert
        if mask_dir == 1:
            # direction small to large
            source = mask_min_vert
            target = mask_max_vert
        return {"source" : source, "target" : target}
    
    def _mask_adj_info(self, mask_vert : int):
        if mask_vert & self._mask_all_verts == 0:
            raise ValueError(f"Invalid mask_vert.")
        succ_arrs = 0
        succ_verts = 0
        pred_arrs = 0
        pred_verts = 0
        for mask_arr in self._mask_arrs.values():
            # check if the vertex is either the source or target of the arrow
            if mask_vert & mask_arr != 0:
                arr_info = self._arr_info_from_mask(mask_arr)
                mask_source = arr_info["source"]
                mask_target = arr_info["target"]
                if mask_source == mask_vert:
                    succ_arrs |= mask_arr
                    succ_verts |= mask_target
                if mask_target == mask_vert:
                    pred_arrs |= mask_arr
                    pred_verts |= mask_source
        pred_info = {"verts" : pred_verts, "arrs" : pred_arrs}
        succ_info = {"verts" : succ_verts, "arrs" : succ_arrs}
        return {"pred_info" : pred_info, "succ_info" : succ_info}
    
    def _mask_pred_info(self, mask_vert : int):
        return self._mask_adj_info(mask_vert)["pred_info"]
    
    def _mask_pred_verts(self, mask_vert : int):
        return self._mask_pred_info(mask_vert)["verts"]
    
    @cached_property
    def _mask_source_verts(self):
        mask_verts = self._mask_verts
        mask_source = 0
        for mask in mask_verts.values():
            if self._mask_pred_verts(mask) == 0:
                mask_source |= mask
        return mask_source
    
    def _mask_pred_arrs(self, mask_vert : int):
        return self._mask_pred_info(mask_vert)["arrs"]
    
    def _mask_succ_info(self, mask_vert : int):
        return self._mask_adj_info(mask_vert)["succ_info"]
    
    def _mask_succ_verts(self, mask_vert : int):
        return self._mask_succ_info(mask_vert)["verts"]
    
    @cached_property
    def _mask_sink_verts(self):
        mask_verts = self._mask_verts
        mask_sink = 0
        for mask in mask_verts.values():
            if self._mask_succ_verts(mask) == 0:
                mask_sink |= mask
        return mask_sink
    
    def _mask_succ_arrs(self, mask_vert : int):
        return self._mask_succ_info(mask_vert)["arrs"]
    
    def _mask_anc_info(self, mask_vert : int):
        mask_anc_verts = mask_vert
        mask_anc_arrs = 0
        mask_test = mask_vert
        mask_verts = self._mask_verts
        while mask_test != 0:
            for mask in mask_verts.values():
                if mask_test & mask != 0:
                    mask_pred = self._mask_pred_info(mask)
                    mask_anc_verts |= mask_pred["verts"]
                    mask_anc_arrs |= mask_pred["arrs"]
                    mask_test = mask_test & (~ mask) | mask_pred["verts"]
        return {"verts" : mask_anc_verts, "arrs" : mask_anc_arrs}
    
    def _mask_des_info(self, mask_vert : int):
        mask_des_verts = mask_vert
        mask_des_arrs = 0
        mask_test = mask_vert
        mask_verts = self._mask_verts
        while mask_test != 0:
            for mask in mask_verts.values():
                if mask_test & mask != 0:
                    mask_succ = self._mask_succ_info(mask)
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
        return sum(self._mask_anc_verts(v) for v in self._mask_verts if v & mask_verts != 0)
    
    def _mask_des_verts_list(self, mask_verts : list[int]):
        return sum(self._mask_des_verts(v) for v in self._mask_verts if v & mask_verts != 0)

    def _mask_anc_verts_comp(self, mask_vert_list : int):
        comp = []
        mask_verts = self._mask_verts
        mask_source = self._mask_source_verts
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
        mask_verts = self._mask_verts
        mask_sink = self._mask_sink_verts
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
        return len(self._mask_des_verts_comp(self._mask_all_verts)) == 1

    @cached_property
    def _mask_anc_closed(self):
        result = []
        mask_verts = self._mask_verts
        for n in range(1,2**len(mask_verts)):
            test_mod = n << len(self.arrows) + 1
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
        mask_verts = self._mask_verts
        for n in range(1,2**len(mask_verts)):
            test_mod = n << len(self.arrows) + 1
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


        

if __name__ == "__main__":

        # bg = BitmaskGraph([Arrow(0, 1, "a"), Arrow(1, 2, "b"), Arrow(2, 3, "c")])
    # mask_verts = bg._mask_verts
    # for vert, mask in mask_verts.items():
    #     print(vert, mask, bin(mask))
    # for arr, mask in bg._mask_arrs.items():
    #     print(arr, bin(mask_verts[arr.source]), bin(mask_verts[arr.target]), mask, bin(mask))

    # print(bg._mask_all_verts & 16)
    # print(mask_verts.values())

    # for i in mask_verts.values():
    #     print(i, bg._mask_adj_info(i))
    #     print(bg._mask_anc_info(i))

    bg1 = BitmaskGraph([Arrow(0,1), Arrow(1,2), Arrow(2,3)])
    vert2 = bg1._mask_all_verts
    mask = bg1._mask_verts
    print(bg1._mask_info)
    # print("Quotient modules")
    # for i in bg1._mask_anc_verts_comp(mask[2] | mask[4]):
    #     print(i, bg1._mask_to_verts(i))
    # print("submodules")
    # for i in bg1._mask_des_verts_comp(mask[2] | mask[4]):
    #     print(i, bg1._mask_to_verts(i))
    print("=== SUBMODS ===")
    des_closed = bg1.des_closed
    for i in des_closed:
        print(i)

    print("=== QUOTIENTS ===")
    anc_closed = bg1.anc_closed
    for i in anc_closed:
        print(i)

    print("=== NEW ===")
    anc_closed_from_des = bg1.anc_closed_from_des
    for i in anc_closed_from_des:
        print(i)
    




