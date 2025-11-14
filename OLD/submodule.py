from networkx.classes import subgraph
import networkx as nx

class Submodule:
    def __init__(self, super_vertex_labels:list, super_arrows:list, sub_vertex_labels:list, inclusion_map:dict, arrow_labels=None):
        self.super_vert = super_vertex_labels
        self.super_arr = super_arrows
        self.sub_vert = sub_vertex_labels
        self.inclusion_map = inclusion_map
        self.sub_arr = self.arrows_of_submodule()
        self.arrow_labels = arrow_labels

    def arrows_of_submodule(self):
        sub_arr = []
        for (x,y) in self.super_arr:
            if x and y in self.sub_vert:
                sub_arr.append((x,y))
        return sub_arr


if __name__ == "__main__":

    S = Submodule(['0','1','2'],[(0,1),(2,0),(2,1),(1,1)], ['0','2'], {0:0,1:2})
    print(S.sub_arr)
