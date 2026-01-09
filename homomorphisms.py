from bitmask_subgraph import Vertex, Arrow
from modulediagram import ModuleDiagram, ModuleSubDiagram
import time
import numpy as np
from numpy.typing import NDArray
import networkx as nx

class Homomorphism:
    def __init__(self,
                 domain: ModuleDiagram,
                 codomain: ModuleDiagram,
                 mapping: dict[Vertex, Vertex]): 
        self.domain = domain
        self.codomain = codomain
        self.mapping = mapping
        self.index_mapping = {self.domain.vertex_to_index[k] : self.codomain.vertex_to_index[v]
                                for k, v in self.mapping.items()}
        self._check_validity()
        
        self.dimension = len(self.mapping)

        

    def _check_validity(self):
        """
        Check homomorphism is valid. Must satisfy:
            - Keys of mapping form a quotient diagram of domain.
            - Values of mapping form a subdiagram of codomain.
            - Quotient and subdiagram are isomorphic.
        """
        for k, v in self.mapping.items():
            if k not in self.domain.vertex_list:
                raise ValueError(f"Invalid mapping. Key {k} not a vertex index of domain.")
            if v not in self.codomain.vertex_list:
                raise ValueError(f"Invalid mapping. Value {v} not a vertex index of codomain.")
        dom_im = ModuleSubDiagram(self.domain, list(self.mapping.keys()))
        if not dom_im.is_quotient():
            raise ValueError(f"The keys of the mapping are not a quotient diagram of the domain.")
        codom_im = ModuleSubDiagram(self.codomain, list(self.mapping.values()))
        if not codom_im.is_submodule():
            raise ValueError(f"The values of the mapping are not a subdiagram of the codomain.")
        quotient = ModuleSubDiagram(self.domain, list(self.mapping.keys()))
        sub = ModuleSubDiagram(self.codomain, list(self.mapping.values()))
        if quotient != sub:
            raise ValueError(f"The image of the homomorphism is not well-defined.")

        self.domain = domain
        self.image = image
        self.codomain = codomain

    