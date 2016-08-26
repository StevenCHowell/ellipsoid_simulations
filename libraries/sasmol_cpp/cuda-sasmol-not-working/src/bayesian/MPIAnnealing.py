# -*- coding: utf-8 -*-


import altar, pyre


class MPIAnnealing(pyre.component, family='altar.bayesian.MPIAnnealing', implements=altar.protocols.Worker):


    samples     = pyre.properties.int()
    parameters  = pyre.properties.int()
    
    worker      = pyre.properties.facility(protocol = altar.protocols.Worker)


    def __init__(self, *, communicator, **kwds):
        super().__init__(**kwds)
        
        self.communicator = communicator
        
        self.ext = self.XMPIAnnealing(
            self.communicator.capsule,
            self.samples,
            self.parameters,
            self.worker.ext,
            )
        
        return


    # private data
    from .ext import XMPIAnnealing # my extension type


# end of file 
