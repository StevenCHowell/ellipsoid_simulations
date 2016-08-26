# -*- coding: utf-8 -*-


import altar, pyre


class Annealer(pyre.component, family='altar.bayesian.Annealer', implements=altar.protocols.Annealer):


    problem     = pyre.properties.facility(protocol = altar.protocols.Problem)
    scheduler   = pyre.properties.facility(protocol = altar.protocols.Scheduler)
    sampler     = pyre.properties.facility(protocol = altar.protocols.Sampler)
    worker      = pyre.properties.facility(protocol = altar.protocols.Worker)

    chains      = pyre.properties.int(default = 2)
    tolerance   = pyre.properties.float(default = 0.005)


    def __init__(self, **kwds):
        super().__init__(**kwds)
        
        self.ext = self.XAnnealer(
            self.problem.ext,
            self.scheduler.ext,
            self.sampler.ext,
            self.worker.ext,
            self.chains,
            self.tolerance,
            )
        
        return


    # interface
    @pyre.export
    def posterior(self, restart = 0, wfreq = 0, fprior = ""):
        self.ext.posterior(restart, wfreq, fprior);
        return


    # private data
    from .ext import XAnnealer # my extension type


# end of file 
