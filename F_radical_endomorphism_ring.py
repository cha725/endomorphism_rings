from A_bitmask_subgraph import Vertex, Arrow
from B_quiver_algebra import MonomialQuiverAlgebra
from D_projective_module_diagram import ProjectiveDiagram, PathModuleDiagram
from E_homomorphisms import EndoRing

class RadicalEndo:
    """
    Represents the endomorphism ring of the radical submodules of an algebra.

    Attributes:
        - quiver_algebra (MonomialQuiverAlgebra): the algebra to take the endomorphism
            ring over and that defines the radical submodules.
    """
    def __init__(self,
                 quiver_algebra: MonomialQuiverAlgebra):
        self.quiver_algebra = quiver_algebra

    def _create_radical_endo(self):
        """
        1. Create the projective modules over the quiver algebra.
        2. Compute their indecomposable summands.
            This will be the module to take the endomorphisms of.
        3. Compute their endomorphism ring.
        """
        projs = [ProjectiveDiagram(self.quiver_algebra, v) for v in self.quiver_algebra.vertices]
        rad_submods = []

        print("\n== Projectives ==")
        for proj in projs:
            poss_submods = proj.radical_submodules
            for p in poss_submods:
                p = PathModuleDiagram(p.vertex_list,p.arrow_list)
                if p not in rad_submods:
                    rad_submods.append(p)

        print(f"\n== Radical Submodules ==")
        for rad_submod in rad_submods:
            print(rad_submod.vertex_labels)
        endo = EndoRing(rad_submods)

        print(f"\n== Endomorphisms ==")
        all_homs = endo.all_homs
        num_homs = sum(len(col) for row in all_homs for col in row)
        print(f"There are {num_homs} homomorphisms in total")
        for i, row in enumerate(all_homs):
            for j, homs in enumerate(row):
                if homs:
                    dom = endo.ind_summands[i]
                    codom = endo.ind_summands[j]
                    print(f"\nHom( {dom.vertex_labels}, {codom.vertex_labels} ) has:")
                    for h in homs:
                        print(f"  Hom({h.mapping_str()})")

        print(f"\n== Quiver Arrows ==")
        for arrow in endo.arrows_of_quiver():
            print(f"{arrow.mapping_str()}")
        return endo
    

###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################



if __name__ == "__main__":

    quiver = MonomialQuiverAlgebra([Arrow(Vertex(0),Vertex(1),"a"), 
                                    Arrow(Vertex(1),Vertex(2),"b")])

    RE = RadicalEndo(quiver)
    RE._create_radical_endo()
