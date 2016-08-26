

IF not IMPLEMENTATION:
    
    cdef extern from "altar/bayesian/Annealer.h" namespace "altar::bayesian":
        
        cdef cppclass Annealer:
            
            Annealer(
                Problem problem,
                COV scheduler,
                Metropolis sampler,
                AnnealingMethod worker,
                size_t chains,
                double tolerance,
                ) except +
            
            int posterior(size_t restart, size_t wfreq, string fprior)


cdef class XAnnealer:
    
    IF not IMPLEMENTATION:
        
        cdef Annealer *thisptr   # hold a C++ instance which we're wrapping
        
        cdef XProblem problem
        cdef XCOV scheduler
        cdef XMetropolis sampler
        cdef XAnnealingMethod worker
    
    ELSE:
        
        def __cinit__(self,
                      XProblem problem,
                      XCOV scheduler,
                      XMetropolis sampler,
                      XAnnealingMethod worker,
                      size_t chains,
                      double tolerance,
                      ):
            
            # retain my children
            self.problem = problem
            self.scheduler = scheduler
            self.sampler = sampler
            self.worker = worker
            
            self.thisptr = new Annealer(
                self.problem.thisptr[0],
                self.scheduler.thisptr[0],
                self.sampler.thisptr[0],
                self.worker.annealingMethod[0],
                chains,
                tolerance,
                )
        
        def __dealloc__(self):
            del self.thisptr
    
    
        def posterior(self, size_t restart = 0, size_t wfreq = 0, fprior = ""):
            return self.thisptr.posterior(restart, wfreq, fprior.encode('UTF-8'))
