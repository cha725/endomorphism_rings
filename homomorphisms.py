from bitmask_subgraph import Vertex, Arrow
from modulediagram import ModuleDiagram, ModuleSubDiagram
import time
import numpy as np
from numpy.typing import NDArray
import networkx as nx

class Homomorphism:
    def __init__(self,
                 domain : ModuleDiagram,
                 image : ModuleDiagram,
                 codomain : ModuleDiagram):
        
        
        self.domain = domain
        self.image = image
        self.codomain = codomain

    