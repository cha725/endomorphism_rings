from modulediagram import ModuleDiagram

class Homomorphism:
    def __init__(self,
                 domain : ModuleDiagram,
                 image : ModuleDiagram,
                 codomain : ModuleDiagram):
        
        
        self.domain = domain
        self.image = image
        self.codomain = codomain

    