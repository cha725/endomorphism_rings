import networkx as nx
import time

from functools import cached_property

from A_bitmask_subgraph import Vertex, Arrow
from C_module_diagram import ModuleDiagram, ModuleSubDiagram

class Homomorphism:
    """
    Represents a homomorphism between two module diagrams by a dictionary
    mapping between the vertices.

    Attributes:
        - domain (ModuleDiagram): The source of the homomorphism
            i.e. the keys in the mapping
        - codomain (ModuleDiagram): The target of the homomorphism
            i.e. the values in the mapping
        - mapping (dict[Vertex, Vertex]): The mapping from vertices in the domain
            to vertices in the codomain.
    """
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
        if not self.domain.is_quotient(list(self.mapping.keys())):
            raise ValueError(f"The keys of the mapping are not a quotient diagram of the domain.")
        if not self.codomain.is_submodule(list(self.mapping.values())):
            raise ValueError(f"The values of the mapping are not a subdiagram of the codomain.")
        quotient = ModuleSubDiagram(self.domain, list(self.mapping.keys()))
        sub = ModuleSubDiagram(self.codomain, list(self.mapping.values()))
        if not quotient.compute_isomorphism(sub):
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
    
    def is_injective(self) -> bool:
        """ 
        Check if the map is injective. 
        i.e. every vertex of the domain is a key in the mapping. 
        """
        return set(self.mapping.keys()) == set(self.domain.vertex_list)
    
    def is_surjective(self) -> bool:
        """ 
        Check if the map is surjective. 
        i.e. every vertex of the codomain is a value in the mapping.
        """
        return set(self.mapping.values()) == set(self.codomain.vertex_list)
    
    def post_compose(self, other: "Homomorphism") -> "Homomorphism | None":
        """
        Post-compose homomorphism with other, -self-> -other->.

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
    
    def pre_compose(self, other: "Homomorphism") -> "Homomorphism | None":
        """
        Pre-compose homomorphism with other, -other-> -self->.

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
    
    def is_isomorphism(self) -> bool:
        """ Check if the map is a bijection between the domain and codomain. """
        if not self.is_injective():
            return False
        if not self.is_surjective():
            return False
        return self.domain == self.codomain
    
    def __eq__(self, other: "Homomorphism"):
        dom = (self.domain == other.domain)
        codom = (self.codomain == other.codomain)
        map = (self.mapping.items() == other.mapping.items())
        return all([dom,codom,map])
    
    def __hash__(self):
        return hash((self.domain, self.codomain, self.mapping))

    def mapping_str(self):
        return "{" + ", ".join(f"{k.label} : {v.label}" for k, v in self.mapping.items()) + "}"

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
                    mapping = q_mod.compute_isomorphism(s_mod)
                    if mapping:
                        try:
                            hom = Homomorphism(self.domain, self.codomain, mapping)
                            homs.append(hom)
                        except:
                            pass
        return homs
    
    def __repr__(self):
        return f"Hom({self.domain.vertex_labels, self.codomain.vertex_labels})"

class EndoRing:
    """
    Represents the endomorphism ring of a direct sum of modules.

    Takes in a list of module diagrams and computes all morphisms between them.
    Finds the indecomposable morphisms and produces the quiver of the endomorphism
    ring.

    Attributes:
        - modules (list[ModuleDiagram]): list of module diagrams whose direct sum to take
        - cut_off (int): maximum number of compositions to take in endomorphism ring.
            If cut_off is low, then many of the endomorphisms could be missing.
    """

    def __init__(self,
                 modules: list[ModuleDiagram],
                 cut_off: int = 50):
        self.modules = modules
        self.cut_off = cut_off
        ind_summands: list[ModuleSubDiagram] = []
        for m in self.modules:
            ind_summands += m.indecomposable_summands
        self.ind_summands: tuple[ModuleSubDiagram,...] = tuple(set(ind_summands))
        self.summand_to_index: dict[ModuleDiagram,int] = {s : idx for idx, s in enumerate(self.ind_summands)}
        self.num_summands: int = len(self.ind_summands)

    def quiver(self):
        """
        Create the quiver of the endomorphism ring as a NetworkX MultiDiGraph.
        The vertices correspond to indecomposable summands of the module M in End(M).
        The arrows are the indecomposable morphisms between summands of M and not the
        isomorphisms which represent the stationary paths.
        """
        G = nx.MultiDiGraph()
        vertices = self.vertices_of_quiver()
        for v in vertices:
            G.add_node(v)
        arrows = self.arrows_of_quiver()
        for arrow in arrows:
            if not arrow.is_isomorphism():
                G.add_edge(arrow.domain,arrow.codomain)
        return G
    
    def draw_quiver(self, node_size=800, node_color='lightblue', edge_color='black'):
        """
        Draw the quiver of the endomorphism ring.
        The vertices correspond to indecomposable summands of the module M in End(M).
        The arrows are the indecomposable morphisms between summands of M.
        """
        G = self.quiver()
        pos = nx.spring_layout(G, seed=42)
        labels = {n: str(n.vertex_list) for n in G.nodes()}
        nx.draw(
            G, pos,
            node_size=node_size,
            node_color=node_color,
            edge_color=edge_color,
            with_labels=True,
            labels=labels
        )
    
    def vertices_of_quiver(self):
        """
        The vertices of the quiver of the endomorphism ring End(M) are precisely the
        indecomposable summands of M.
        """
        return self.ind_summands
    
    def arrows_of_quiver(self) -> list[Homomorphism]:
        """
        The arrows of the quiver of the endomorphism ring End(M) are precisely the 
        indecomposable morphisms between indecomposable summands of M.
        i.e. the morphisms f such that if f=gh then either g or h is idempotent.
        """
        indecomp_homs = self.find_indecomposable_morphisms()
        arrows = []
        for row in indecomp_homs:
            for homs in row:
                for hom in homs:
                    arrows.append(hom)
        return arrows

    def find_indecomposable_morphisms(self) -> list[list[list[Homomorphism]]]:
        """
        Identify all homomorphisms between indecomposable summands that do not 
        factor non-trivially through another homomorphism.

        i.e. find all the morphisms f such that if f = gh then either g or h is
            idempotent.
        
        Returns:
            - list[list[list[Homomorphism]]]:
                A nested list (like the input for a numpy matrix) whose (i,j)-entry
                is the list of indecomposable morphisms from the ith indecomposable
                summand to the jth.
        
        TODO: Is there a better way to find these morphisms that does not require
        computing all endomorphisms first?
        """
        candidate_homs = self.all_homs
        indecomposables = [[list(h) for h in row] for row in self.all_homs]
        for row in range(self.num_summands):
            for col in range(self.num_summands):
                for hom in candidate_homs[row][col]:
                    if not self.is_indecomposable_hom(hom):
                        indecomposables[row][col].remove(hom)
        return indecomposables
    
    def is_indecomposable_hom(self, hom: Homomorphism) -> bool:
        """ Check whether given hom is indecomposable. """
        dom = hom.domain
        codom = hom.codomain
        try:
            composed_homs = self.composed_homs(dom, codom)
        except Exception as e:
            raise e
        for composed in composed_homs[1:]:
            if hom in composed:
                return False
            if hom.is_isomorphism():
                return False
        return True
    
    def composed_homs(self, domain: ModuleDiagram, codomain: ModuleDiagram) -> list[list[Homomorphism]]:
        """ 
        Return list of lists of homomorphisms where the ith entry of the list are the
        homomorphisms domain to codomain that are compositions of i homs.
        """
        try:
            dom_idx = self.summand_to_index[domain]
        except:
            raise ValueError(f"Domain {domain} is not an indecomposable summand of M.")
        try:
            codom_idx = self.summand_to_index[codomain]
        except:
            raise ValueError(f"Codomain {codomain} is not an indecomposable summand of M.")
        composed_homs : list[list[Homomorphism]] = []
        all_composed_homs = self.all_composed_homs()
        for n_composed_homs in all_composed_homs:
            composed_homs.append(n_composed_homs[dom_idx][codom_idx])
        return composed_homs
    
    def all_composed_homs(self, num_compositions: int | None = None) -> list[list[list[list[Homomorphism]]]]:
        """
        Return list of compositions of all homs in starting homs up to n compositions.
        """
        num_compositions = num_compositions or self.cut_off
        previous_homs = self.all_homs
        comp_homs = [self.all_homs]
        for _ in range(2, num_compositions + 1):
            new_homs = self._compose_homs(previous_homs, self.all_homs)
            if new_homs == [[[] for _ in range(self.num_summands)] for _ in range(self.num_summands)]:
                break
            comp_homs.append(new_homs)
            previous_homs = new_homs
        return comp_homs
    
    def _compose_homs(self, previous_homs: list[list[list[Homomorphism]]], new_homs: list[list[list[Homomorphism]]]) -> list[list[list[Homomorphism]]]:
        """ Return composition of previous_homs then new_homs. """
        comp_homs: list[list[list[Homomorphism]]] = [[[] for col in range(self.num_summands)] for row in range(self.num_summands)]
        for i in range(self.num_summands):
            for j in range(self.num_summands):
                seen = []
                for k in range(self.num_summands):
                    # find compositions i -f-> k -g-> j
                    for f in previous_homs[i][k]:
                        if f.is_isomorphism():
                            continue
                        for g in new_homs[k][j]:
                            if g.is_isomorphism():
                                continue
                            try:
                                comp_hom = g.pre_compose(f)
                                if comp_hom in seen:
                                    continue
                                if comp_hom:
                                    seen.append(comp_hom)
                                    comp_homs[i][j].append(comp_hom)
                            except:
                                continue
        return comp_homs
                    
    @cached_property
    def all_homs(self) -> list[list[list[Homomorphism]]]:
        """ 
        Compute all homomorphisms between indecomposable summands.

        Returns:
            - list[list[list[Homomorphism]]]:
                A nested list (like the input for a numpy matrix) whose (i,j)-entry 
                is the list of homomorphisms from the ith indecomposable summand to 
                the jth.
        """
        return [[HomomorphismGroup(m, n).homs for n in self.ind_summands] for m in self.ind_summands]


    def generate_homs_from_ind(self, num_compositions: int | None = None) -> list[list[list[Homomorphism]]]:
        num_compositions = num_compositions or self.cut_off
        all_previous_homs = [self.find_indecomposable_morphisms()]
        ind_homs = self.find_indecomposable_morphisms()
        for _ in range(2, num_compositions + 1):
            previous_homs = all_previous_homs[-1]
            new_homs = self._compose_homs(previous_homs, ind_homs)
            if new_homs == [[[] for col in range(self.num_summands)] for row in range(self.num_summands)]:
                break
            all_previous_homs.append(new_homs)
        all_homs = [[[] for _ in range(self.num_summands)] for _ in range(self.num_summands)]
        for dom in range(self.num_summands):
            for codom in range(self.num_summands):
                homs_dom_to_codom: list[Homomorphism] = []
                for matrix in all_previous_homs:
                    for hom in matrix[dom][codom]:
                        if hom not in homs_dom_to_codom:
                            homs_dom_to_codom.append(hom)
                all_homs[dom][codom] = homs_dom_to_codom
        return all_homs
    
    def __repr__(self):
        return f"EndoRing(modules={self.modules})"
    





        
    
class ExamplesHomomorphisms:
    """
    Collection of example module diagrams and homomorphisms.
    Ensures all examples are valid homomorphisms according to the rules in Homomorphism.
    """
    def __init__(self):
        self.examples = {}

    def add(self, 
            name: str, 
            domain: ModuleDiagram,
            codomain: ModuleDiagram,
            mapping: dict[Vertex, Vertex]):
        self.examples[name] = (domain, codomain, mapping)


    def add_simple_chain(self, name: str):
        """
        Add a canonical example:
        Domain:   0 -> 1 -> 2
        Codomain: 0 -> 1 -> 2 -> 3
        Mapping:  0->1, 1->2, 2->3
        """
        dom = ModuleDiagram(arrow_list=[Arrow(Vertex(0),Vertex(1)),
                                    Arrow(Vertex(1),Vertex(2))]
                            )
        cod = ModuleDiagram(arrow_list=[Arrow(Vertex(0),Vertex(1)),
                                    Arrow(Vertex(1),Vertex(2)),
                                    Arrow(Vertex(2),Vertex(3))]
                            )
        mapping = {Vertex(0):Vertex(1), Vertex(1):Vertex(2), Vertex(2):Vertex(3)}
        self.add(name, dom, cod, mapping)

    def add_invalid_simple_chain(self, name: str):
        """
        Add a canonical example:
        Domain:   0 -> 1 -> 2
        Codomain: 0 -> 1 -> 2 -> 3
        Mapping:  0->0, 1->1, 2->2
        """
        dom = ModuleDiagram(
            arrow_list = [
                Arrow(Vertex(0),Vertex(1)),
                Arrow(Vertex(1),Vertex(2))
            ]
        )
        cod = ModuleDiagram(
            arrow_list = [
                Arrow(Vertex(0),Vertex(1)),
                Arrow(Vertex(1),Vertex(2)),
                Arrow(Vertex(2),Vertex(3))
            ]
        )
        mapping = {Vertex(0):Vertex(0), Vertex(1):Vertex(1), Vertex(2):Vertex(2)}
        self.add(name, dom, cod, mapping)

    def add_projection(self, name: str):
        """
        Domain:   0 -> 1 -> 2
        Codomain: 0 -> 1
        Mapping:  0->0, 1->1   (drop vertex 2)
        """
        dom = ModuleDiagram(
            arrow_list = [
                Arrow(Vertex(0),Vertex(1)),
                Arrow(Vertex(1),Vertex(2))
            ]
        )
        cod = ModuleDiagram(
            arrow_list = [
                Arrow(Vertex(0),Vertex(1))
            ]
        )
        mapping = {Vertex(0):Vertex(0), Vertex(1):Vertex(1)}
        self.add(name, dom, cod, mapping)

    def add_inclusion(self, name: str):
        """
        Domain:   0 -> 1
        Codomain: 0 -> 1 -> 2
        Mapping:  0->1, 1->2   (include as a subdiagram)
        """
        dom = ModuleDiagram(
            arrow_list = [
                Arrow(Vertex(0),Vertex(1))
            ]
        )
        cod = ModuleDiagram(
            arrow_list = [
                Arrow(Vertex(0),Vertex(1)),
                Arrow(Vertex(1),Vertex(2))
            ]
        )
        mapping = {Vertex(0):Vertex(1), Vertex(1):Vertex(2)}
        self.add(name, dom, cod, mapping)

    def add_invalid_inclusion(self, name: str):
        """
        Domain:   0 -> 1
        Codomain: 0 -> 1 -> 2
        Mapping:  0->0, 1->1   (include as a subdiagram)
        """
        dom = ModuleDiagram(
            arrow_list = [
                Arrow(Vertex(0),Vertex(1))
            ]
        )
        cod = ModuleDiagram(
            arrow_list = [
                Arrow(Vertex(0),Vertex(1)),
                Arrow(Vertex(1),Vertex(2))
            ]
        )
        mapping = {Vertex(0):Vertex(0), Vertex(1):Vertex(1)}
        self.add(name, dom, cod, mapping)

    def run(self, verbose=True):
        results = []

        for name, (dom, cod, mapping) in self.examples.items():
            print(f"\n=== Homomorphism Example: {name} ===")

            print(f"dom = {dom}")
            print(f"cod = {cod}")
            print(f"\nmapping = {mapping}.")

            try:
                start = time.time()
                h = Homomorphism(dom, cod, mapping)
                init_time = time.time() - start
                print(f"\nConstruction time: {init_time:.5f}s")

            except Exception as e:
                print(f"\nInvalid homomorphism: {e}")
                continue

            # Image
            start = time.time()
            img = h.image()
            img_time = time.time() - start
            print(f"Image computation time: {img_time:.5f}s")

            # Test composition with identity (when possible)
            if all(v in mapping.values() for v in mapping.keys()):
                # Build identity-like map on codomain
                id_map = {v: v for v in mapping.values()}
                try:
                    g = Homomorphism(cod, cod, id_map)
                    comp = h.post_compose(g)
                    print("Composition with identity succeeded.")
                except Exception:
                    print("Composition failed (identity check).")

            if verbose:
                print("\nDomain diagram vertices:", dom.vertex_list)
                print("Codomain diagram vertices:", cod.vertex_list)
                print("Mapping:", mapping)
                print("Image vertices:", img.vertex_list)

            results.append((init_time, img_time))

        return results
    
class ExamplesHomGroups:
    """
    Collection of example (domain, codomain) pairs for computing Hom(domain, codomain).
    Uses HomomorphismGroup to enumerate all homomorphisms.
    """

    def __init__(self):
        self.examples = {}

    def add(self, name: str,
                  domain: ModuleDiagram,
                  codomain: ModuleDiagram):
        """
        Add a named hom-group example consisting of (domain, codomain).
        """
        self.examples[name] = (domain, codomain)

    def add_simple_chain_pair(self, name: str):
        """
        Domain:   0 -> 1 -> 2
        Codomain: 0 -> 1 -> 2 -> 3
        """
        dom = ModuleDiagram(
            arrow_list = [
                Arrow(Vertex(0),Vertex(1)),
                Arrow(Vertex(1),Vertex(2))
            ]
        )
        cod = ModuleDiagram(
            arrow_list = [
                Arrow(Vertex(0),Vertex(1)),
                Arrow(Vertex(1),Vertex(2)),
                Arrow(Vertex(2),Vertex(3))
            ]
        )
        self.add(name, dom, cod)

    def add_identity_pair(self, name: str):
        """
        Domain = Codomain = 0 -> 1 -> 2
        """
        M = ModuleDiagram(
            arrow_list = [
                Arrow(Vertex(0),Vertex(1)),
                Arrow(Vertex(1),Vertex(2))
            ]
        )
        self.add(name, M, M)

    def add_projection_pair(self, name: str):
        """
        Domain:   0 -> 1 -> 2
        Codomain: 0 -> 1
        """
        dom = ModuleDiagram(
            arrow_list = [
                Arrow(Vertex(0),Vertex(1)),
                Arrow(Vertex(1),Vertex(2))
            ]
        )
        cod = ModuleDiagram(
            arrow_list = [
                Arrow(Vertex(0),Vertex(1))
            ]
        )
        self.add(name, dom, cod)

    def add_inclusion_pair(self, name: str):
        """
        Domain:   0 -> 1
        Codomain: 0 -> 1 -> 2
        """
        dom = ModuleDiagram(
            arrow_list = [
                Arrow(Vertex(0),Vertex(1))
            ]
        )
        cod = ModuleDiagram(
            arrow_list = [
                Arrow(Vertex(0),Vertex(1)),
                Arrow(Vertex(1),Vertex(2))
            ]
        )
        self.add(name, dom, cod)

    def run(self, verbose=True):
        """
        Compute Hom(domain, codomain) for every stored example.
        Returns a dict: name → (number_of_homs, computation_time).
        """
        results = {}

        for name, (dom, cod) in self.examples.items():
            print(f"\n=== HomGroup Example: {name} ===")

            start = time.time()
            HG = HomomorphismGroup(dom, cod)
            homs = HG.compute_homs()
            elapsed = time.time() - start

            print(f"Computed {len(homs)} homomorphisms in {elapsed:.5f}s.")

            if verbose:
                print(f"Domain = {dom}")
                print(f"Codomain = {cod}")
                for h in homs:
                    print(h)

            results[name] = (len(homs), elapsed)

        return results

class ExamplesEndoRing:
    """
    Collection of example module diagrams for testing EndoRing(M).
    Runs the construction and optional indecomposable–morphism search.
    """
    def __init__(self):
        self.examples = {}

    def add(self, name: str, M: list[ModuleDiagram]):
        print("add", M)
        self.examples[name] = M

        

    def add_chain(self, name: str):
        """
        Example: 0 -> 1 -> 2 -> 3
        """
        M = ModuleDiagram(
            arrow_list = [
                Arrow(Vertex(0),Vertex(1)),
                Arrow(Vertex(1),Vertex(2)),
                Arrow(Vertex(2),Vertex(3))
            ]
        )
        self.add(name, [M])

    def add_short_chain(self, name: str):
        """
        Example: 0 -> 1 -> 2
        """
        M = ModuleDiagram(
            arrow_list = [
                Arrow(Vertex(0),Vertex(1)),
                Arrow(Vertex(1),Vertex(2))
            ]
        )
        self.add(name, [M])

    def add_double_fork(self, name: str):
        """
        Example: 0 -> 1 <- 2 -> 3
        """
        M = ModuleDiagram(
            arrow_list = [
                Arrow(Vertex(0),Vertex(1)),
                Arrow(Vertex(2),Vertex(1)),
                Arrow(Vertex(2),Vertex(3))
            ]
        )
        self.add(name, [M])

    def add_custom(self, 
                   name: str,
                   arrow_list: list[Arrow]):
        """
        Add user-defined module.
        """
        M = ModuleDiagram(
            arrow_list=arrow_list
        )
        self.add(name, [M])

    def run(self, verbose=True, compute_indecomposables=False, cut_off=50):
        """
        For each example M, construct EndoRing(M) and optionally
        compute indecomposable morphisms.
        """
        results = {}

        for name, M in self.examples.items():
            print(f"\n=== EndoRing Example: {name} ===")
            print(f"Module: {M}")

            #try:
            start = time.time()
            ER = EndoRing(M, cut_off=cut_off)
            elapsed = time.time() - start
            print(f"EndoRing constructed in {elapsed:.5f}s")

            if verbose:
                print(f"Indecomposable summands: {ER.ind_summands}")
                print(f"No. summands: {ER.num_summands}")
                print(f"\nElements:")
                for i, row in enumerate(ER.all_homs):
                        for j, homs in enumerate(row):
                            print(f"  Hom({i},{j}) has:")
                            for h in homs:
                                print(h)

            indecomp = None
            if compute_indecomposables:
                print("\nComputing indecomposable morphisms...")
                s = time.time()
                indecomp = ER._find_indecomposable_morphisms()
                
                
                t = time.time() - s
                print(f"\nFound indecomposable maps in {t:.5f}s")
                if verbose:
                    for i, row in enumerate(indecomp):
                        for j, homs in enumerate(row):
                            print(f"  Hom({i},{j}) has:")
                            for h in homs:
                                print(h)

            print("\ncompute homs as comps of inds")
            new_gens = ER.find_ind_homs_as_comp(ER.cut_off)
            for h in new_gens:
                print(h)
            print("\n check generate:", len(ER.all_homs), len(new_gens))

            quiver = ER.quiver()
            print(quiver)

            results[name] = {
                "construction_time": elapsed,
                "indecomposables": indecomp
            }

            #except Exception as e:
              #  print(f"EndoRing failed: {e}")
        return results





if __name__ == "__main__":

    # # ---------------------------------------------------------
    # # Examples of single homomorphisms
    # # ---------------------------------------------------------

    # examples = ExamplesHomomorphisms()
    # examples.add_inclusion("inc")
    # examples.add_projection("proj")
    # examples.add_simple_chain("a2")
    # examples.run(verbose=True)

    # # ---------------------------------------------------------
    # # Examples of hom groups
    # # ---------------------------------------------------------

    # print("\n=== Hom Groups ===")

    # hom_groups = ExamplesHomGroups()
    # hom_groups.add_simple_chain_pair("chain_pair")
    # hom_groups.add_identity_pair("id_pair")
    # hom_groups.add_projection_pair("proj_pair")
    # hom_groups.add_inclusion_pair("inc_pair")

    # results = hom_groups.run(verbose=True)

    # print("\nSummary:")
    # for name, (n, t) in results.items():
    #     print(f"{name}: {n} homs in {t:.5f}s")


    # print("\n=== EndoRing Examples ===")

    examples = ExamplesEndoRing()
    examples.add_chain("chain4")
    examples.add_short_chain("chain3")
    examples.add_double_fork("fork")
    examples.add_custom(
        "custom",
        arrow_list=[Arrow(Vertex(0),Vertex(1)), Arrow(Vertex(1),Vertex(2))]
    )

    results = examples.run(verbose=True,
                           compute_indecomposables=True,
                           cut_off=40)

    print("\nSummary:")
    for name, data in results.items():
        print(f"{name}: EndoRing built in {data['construction_time']:.5f}s")

    E = EndoRing([ModuleDiagram(arrow_list=[Arrow(Vertex(0),Vertex(1)),Arrow(Vertex(1),Vertex(2)),Arrow(Vertex(2),Vertex(3))]),
             ModuleDiagram(arrow_list=[Arrow(Vertex(0),Vertex(1)),Arrow(Vertex(1),Vertex(2))]),
             ModuleDiagram(arrow_list=[Arrow(Vertex(0),Vertex(1))])])
    print(E.ind_summands)
    indecomp = E._find_indecomposable_morphisms()
    for i, row in enumerate(indecomp):
        for j, homs in enumerate(row):
            print(f"  Hom({i},{j}) has:")
            for h in homs:
                print(h)



