from networkx.classes import subgraph
import networkx as nx

class QuotientModule:
    def __init__(self, super_vertices:list, super_arrows:list, quot_vertices:list, projection_map:dict):
        self.super_vert = super_vertices
        self.super_arr = super_arrows
        self.quot_vert = quot_vertices
        self.projection_map = projection_map

if __name__ == "__main__":

    print('hi')

