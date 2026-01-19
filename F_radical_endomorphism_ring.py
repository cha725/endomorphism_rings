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
        for proj in projs:
            poss_submods = proj.radical_submodules
            for p in poss_submods:
                if p not in rad_submods:
                    rad_submods.append(p)
        endo = EndoRing(rad_submods)
        endo.draw_quiver()
        plt.show()
        indecomp = endo._find_indecomposable_morphisms()
        for i, row in enumerate(indecomp):
            for j, homs in enumerate(row):
                print(f"  Hom({i},{j}) has:")
                for h in homs:
                    print(h)
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

    quiver = MonomialQuiverAlgebra([Arrow(0,1,"a"), 
                                    Arrow(1,2,"b")])

    RE = RadicalEndo(quiver)
    print(RE._create_radical_endo())


    # projs = [ModuleDiagram(arrow_list=[Arrow(Vertex(0),Vertex(1)),Arrow(Vertex(1),Vertex(2)),Arrow(Vertex(2),Vertex(3))]),
    #          ModuleDiagram(arrow_list=[Arrow(Vertex(0),Vertex(1)),Arrow(Vertex(1),Vertex(2))]),
    #          ModuleDiagram(arrow_list=[Arrow(Vertex(0),Vertex(1))])]
    # rad_submods = []
    # for proj in projs:
    #     poss_submods = proj.radical_submodules
    #     for p in poss_submods:
    #         if p not in rad_submods:
    #             rad_submods.append(p)
    # endo = EndoRing(rad_submods)
    # print("created endo")

    # e_quiver = endo.quiver()
    # print(e_quiver)

    # endo.draw_quiver()
    # plt.show()

    # print([s.vertex_list for s in endo.ind_summands])
    # indecomp = endo._find_indecomposable_morphisms()
    # for i, row in enumerate(indecomp):
    #     for j, homs in enumerate(row):
    #         print(f"  Hom({i},{j}) has:")
    #         for h in homs:
    #             print(h)