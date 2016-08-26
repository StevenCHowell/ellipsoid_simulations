# -*- coding: utf-8 -*-


import altar, pyre


class Metropolis(pyre.component, family='altar.bayesian.Metropolis', implements=altar.protocols.Sampler):


    rng = pyre.properties.facility(protocol = altar.protocols.RNG)
    rng.doc = 'the random number generator'
    
    steps = pyre.properties.int(default = 20)
    steps.doc = 'the length of my Markov chains'
    
    scaling = pyre.properties.float(default = .1)
    scaling.doc = 'factor by which the covariance matrix is scaled'
    
    acceptanceWeight = pyre.properties.float(default = 8)
    acceptanceWeight.doc = 'the relative weight of accepted samples'
    
    rejectionWeight = pyre.properties.float(default = 1)
    rejectionWeight.doc = 'the relative weight of rejected samples'


    def __init__(self, **kwds):
        super().__init__(**kwds)
        
        self.ext = self.XMetropolis(
            self.rng.rng.rng,
            self.steps,
            self.scaling,
            self.acceptanceWeight,
            self.rejectionWeight,
            )
        
        return


    # private data
    from .ext import XMetropolis # my extension type


# end of file 
