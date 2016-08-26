# -*- coding: utf-8 -*-


import altar, pyre


class COV(pyre.component, family='altar.bayesian.COV', implements=altar.protocols.Scheduler):

    rng = pyre.properties.facility(protocol = altar.protocols.RNG)

    tolerance       = pyre.properties.float(default = .001)
    maxIterations   = pyre.properties.int(default = 1000)
    target          = pyre.properties.float(default = 1.0)


    def __init__(self, **kwds):
        super().__init__(**kwds)
        
        self.ext = self.XCOV(
            self.rng.rng.rng,
            self.tolerance,
            self.maxIterations,
            self.target,
            )
        
        return


    # private data
    from .ext import XCOV # my extension type


# end of file 
