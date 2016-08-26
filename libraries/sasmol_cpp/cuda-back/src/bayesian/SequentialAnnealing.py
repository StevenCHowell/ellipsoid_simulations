# -*- coding: utf-8 -*-


import altar, pyre


class SequentialAnnealing(pyre.component, family='altar.bayesian.SequentialAnnealing', implements=altar.protocols.Worker):


    samples     = pyre.properties.int()
    parameters  = pyre.properties.int()


    def __init__(self, rank=0, **kwds):
        super().__init__(**kwds)
        
        self.rank = rank
        
        self.ext = self.XSequentialAnnealing(
            self.samples,
            self.parameters,
            self.rank,
            )
        
        return


    # private data
    from .ext import XSequentialAnnealing # my extension type


# end of file 
