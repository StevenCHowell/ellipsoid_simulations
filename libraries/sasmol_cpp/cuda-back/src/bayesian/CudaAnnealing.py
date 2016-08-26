# -*- coding: utf-8 -*-


import altar, pyre


class CudaAnnealing(pyre.component, family='altar.bayesian.CudaAnnealing', implements=altar.protocols.Worker):


    samples     = pyre.properties.int()
    parameters  = pyre.properties.int()


    def __init__(self, *, ngpu, rank=0, **kwds):
        super().__init__(**kwds)

        self.ngpu = ngpu
        self.rank = rank
        
        self.ext = self.XCudaAnnealing(
            self.samples,
            self.parameters,
            self.ngpu,
            self.rank,
            )
        
        return


    # private data
    from .ext import XCudaAnnealing # my extension type


# end of file 
