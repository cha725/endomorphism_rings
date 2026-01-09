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

    def apply(self, vertex: Vertex) -> Vertex:
        """
        Apply the homomorphism to a vertex of the domain module.

        Args:
            vertex (Vertex): A vertex of the domain module.

        Returns:
            int: The image of the element in the codomain module.
        """
        if vertex not in self.mapping:
            raise ValueError(f"Invalid element. Vertex {vertex} not in domain mapping.")
        return self.mapping[vertex]
    
    def image(self) -> ModuleSubDiagram:
        """
        Compute the image of the homomorphism.

        Returns:
            ModuleSubDiagram: subdiagram of the codomain.
        """
        return ModuleSubDiagram(self.codomain, list(self.mapping.values()))
    
    def is_inclusion(self) -> bool:
        return self.image().vertex_labels == self.domain.vertex_labels
    
    def is_projection(self) -> bool:
        return self.image().vertex_labels == self.codomain.vertex_labels
    
        self.domain = domain
        self.image = image
        self.codomain = codomain

    