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
    
    def post_compose(self, other: "Homomorphism") -> Homomorphism | None:
        """
        Post-compose homomorphism with other.

        Args:
            other (Homomorphism): Another homomorphism.

        Returns:
            Homomorphism: The composition of the two mappings in order: self then other.
            None: If the composition is the zero morphism.
        """
        comp_map = {k : u   for k, v in self.mapping.items()
                            for l, u in other.mapping.items()
                            if v == l}
        if comp_map:
            try:
                return Homomorphism(self.domain, other.codomain, comp_map)
            except:
                return None
    
    def pre_compose(self, other: "Homomorphism") -> Homomorphism | None:
        """
        Pre-compose homomorphism with other.

        Args:
            other (Homomorphism): Another homomorphism.

        Returns:
            Homomorphism: The composition of the two mappings in order: other then self.
            None: If the composition is the zero morphism.
        """
        comp_map = {l : v
                    for l, u in other.mapping.items()  
                    for k, v in self.mapping.items()
                    if u == k}
        if comp_map:
            try:
                return Homomorphism(other.domain, self.codomain, comp_map)
            except:
                return None
    
    def post_compose_with_homs(self, post_homs: list[Homomorphism]):
        return [self.post_compose(h) for h in post_homs if self.post_compose(h)]
    
    def pre_compose_with_homs(self, post_homs: list[Homomorphism]):
        return [self.pre_compose(h) for h in post_homs if self.pre_compose(h)]
    
    def is_identity(self):
        if len(self.mapping.keys()) == len(self.domain.vertex_list):
            if len(self.mapping.values()) == len(self.codomain.vertex_list):
                return all(self.mapping[v] == v for v in self.domain.vertex_list)
        return False
    def __eq__(self, other: "Homomorphism"):
        dom = self.domain == other.domain
        codom = self.codomain == other.codomain
        map = self.mapping == other.mapping
        return all([dom,codom,map])
    
    __hash__ = object.__hash__

    def __repr__(self) -> str:
        return f"Hom({self.domain.vertex_labels} -> {self.codomain.vertex_labels} with mapping={dict({k.label : v.label for k,v in self.mapping.items()})})"    

class HomomorphismGroup:
    def __init__(self,
                 domain: ModuleDiagram,
                 codomain: ModuleDiagram):
        self.domain = domain
        self.codomain = codomain
        self.homs = self.compute_homs()

    def compute_homs(self) -> list[Homomorphism]:
        """
        Returns list of homomorphisms from domain to codomain.
        """ 
        homs = []
        quotient_im = self.domain.all_quotients
        sub_im = self.codomain.all_submodules
        max_dim = min(self.domain.num_vertices, self.codomain.num_vertices)
        for dim in range(max_dim + 1):
            poss_quotients = quotient_im[dim]
            poss_subs = sub_im[dim]
            for quotient in poss_quotients:
                for sub in poss_subs:
                    mapping = {}
                    q_mod = ModuleSubDiagram(self.domain, quotient)
                    s_mod = ModuleSubDiagram(self.codomain, sub)
                    if q_mod == s_mod:
                        mapping = q_mod.compute_isomorphism(s_mod)
                        if mapping:
                            homs.append(Homomorphism(self.domain, self.codomain, mapping))
        return homs
    
    def compose_with(self, other: "HomomorphismGroup") -> list[Homomorphism]:
        """
        Return composition of all homs in starting homs up to n compositions.
        """
        if self.codomain != other.domain:
            return []
        new_homs = []
        pre_homs = self.homs
        post_homs = other.homs
        for h in pre_homs:
            new_homs += h.post_compose_with_homs(post_homs)
        return new_homs
    
    def __repr__(self):
        return f"Hom({self.domain.vertex_labels, self.codomain.vertex_labels})"

class EndoRing:

    def __init__(self,
                 modules: list[ModuleDiagram] | ModuleDiagram,
                 cut_off: int = 50):
        self.modules = modules
        self.cut_off = cut_off
        self.ind_summands = []
        if isinstance(modules, list):
            for m in modules:
                self.ind_summands += m.indecomposable_summands()
        else:
            self.ind_summands += modules.indecomposable_summands()
        self.ind_summands = tuple(set(self.ind_summands))
        self.num_summands = len(self.ind_summands)

    def _find_indecomposable_morphisms(self) -> list[list[list[Homomorphism]]]:
        """
        Return list of endomorphisms that do not factor non-trivially through
        another endomorphism. 
        i.e. the endomorphisms f such that if f = gh then either g or h is idempotent.
        """
        candidate_homs = self.all_homs
        indecomposables = [[list(h) for h in row] for row in candidate_homs]

        for i in range(self.num_summands):
            for j in range(self.num_summands):
                new_list = []
                for f in indecomposables[i][j]:
                    reducible = False
                    for k in range(self.num_summands):
                        for g in candidate_homs[i][k]:
                            if g.is_identity():
                                continue
                            for h in candidate_homs[k][j]:
                                if h.is_identity():
                                    continue
                                comp = g.post_compose(h)
                                if comp is not None and comp.hom_signature() == f.hom_signature():
                                    reducible = True
                                    break
                            if reducible:
                                break
                        if reducible:
                            break
                    if not reducible:
                        new_list.append(f)
                indecomposables[i][j] = new_list
        return indecomposables


    def _compose_homs(self, n):
        """
        Return list of compositions of all homs in starting homs up to n compositions.
        """
        num_comps = 0
        starting_homs = self.all_homs.copy()
        composed_homs = [starting_homs.copy()]
        current_homs = self.all_homs.copy()
        while num_comps < n and composed_homs[num_comps] != [[[]]]:
            num_comps += 1
            new_homs: list[list[list[Homomorphism]]] = []
            for i in range(self.num_summands):
                new_homs_from_i: list[list[Homomorphism]] = []
                for j in range(self.num_summands):
                    new_homs_i_to_j: list[Homomorphism] = []
                    for k in range(self.num_summands):
                        pre_homs = current_homs[i][k]
                        post_homs = starting_homs[k][j]
                        for h in pre_homs:
                            if h.is_identity():
                                continue
                            for g in post_homs:
                                if g.is_identity():
                                    continue
                                try:
                                    new_hom = h.pre_compose(g)
                                    if new_hom and new_hom not in new_homs_i_to_j:
                                        new_homs_i_to_j.append(new_hom)
                                except:
                                    continue
                    new_homs_from_i.append(new_homs_i_to_j)
                new_homs.append(new_homs_from_i)
            starting_homs = current_homs
            current_homs = new_homs
            composed_homs.append(current_homs)
        return composed_homs
                    

    def generate_homs_from_ind(self, n: int):
        ind_homs = self._find_indecomposable_morphisms()
        num_comps = 0
        starting_homs = ind_homs
        composed_homs = [starting_homs]
        current_homs = ind_homs
        while num_comps < n and composed_homs[num_comps] != [[[]]]:
            num_comps += 1
            new_homs: list[list[list[Homomorphism]]] = []
            for i in range(self.num_summands):
                new_homs_from_i: list[list[Homomorphism]] = []
                for j in range(self.num_summands):
                    new_homs_i_to_j: list[Homomorphism] = []
                    for k in range(self.num_summands):
                        pre_homs = current_homs[i][k]
                        post_homs = starting_homs[k][j]
                        for h in pre_homs:
                            if h.is_identity():
                                continue
                            for g in post_homs:
                                if g.is_identity():
                                    continue
                                new_hom = h.pre_compose(g)
                                if new_hom and new_hom not in new_homs_i_to_j:
                                    new_homs_i_to_j.append(new_hom)
                    new_homs_from_i.append(new_homs_i_to_j)
                new_homs.append(new_homs_from_i)
            starting_homs = current_homs
            current_homs = new_homs
            composed_homs.append(current_homs)
        return composed_homs
    
    def quiver(self):
        G = nx.DiGraph()
        for m in self.ind_summands:
            G.add_node(m)
        M = self._find_indecomposable_morphisms()
        for row in range(self.num_summands):
            for col in range(self.num_summands):
                for hom in M[row][col]:
                    G.add_edge(hom.domain,hom.codomain)
        return G
    
    def draw_quiver(self, node_size=800, node_color='lightblue', edge_color='black'):
        """
        Draw the quiver of the endomorphism ring.
        The nodes correspond to indecomposable summands of the module M in End(M).
        The arrows are the indecomposable morphisms between summands of M.
        """
        G = self.quiver()
        pos = nx.spring_layout(G, seed=42)

        # Label each node by its vertex list
        labels = {n: str(n.vertex_list) for n in G.nodes()}

        nx.draw(
            G, pos,
            node_size=node_size,
            node_color=node_color,
            edge_color=edge_color,
            with_labels=True,
            labels=labels
        )
        


    