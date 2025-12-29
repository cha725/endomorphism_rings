from quiver_algebra import MonomialQuiverAlgebra
from projective_mod_diagram import ProjectiveDiagram

class RadicalEndo:
    def __init__(self,
                 quiver_algebra: MonomialQuiverAlgebra):
        self.quiver_algebra = quiver_algebra

    def _create_radical_endo(self):
        projs = [ProjectiveDiagram(self.quiver_algebra, v) for v in self.quiver_algebra.vertices]
        rad_submods = [proj.radical_subdiagram() for proj in projs]
        endo = EndomorphismRing(rad_submods)
        return endo

