

IF not IMPLEMENTATION:
    
    cdef extern from "altar/bayesian/cudaMetropolis.h" namespace "altar::bayesian::cuda":
        
        cdef cppclass cudaMetropolis(Metropolis):
            
            cudaMetropolis(
                rng_t *,
                size_t steps,
                double scaling,
                double acceptanceWeight,
                double rejectionWeight,
                ) except +



cdef class XCudaMetropolis(XMetropolis):
    
    IF IMPLEMENTATION:
        
        def __cinit__(self,
                      rngCapsule,
                      size_t steps,
                      double scaling,
                      double acceptanceWeight,
                      double rejectionWeight,
                      ):
            # retain my children
            self.rngCapsule = rngCapsule
            cdef rng_t *rng
            rng = <rng_t *>PyCapsule_GetPointer(self.rngCapsule, rng_capsule_name)
            # XXX: this leaks the Metropolis created in the base class
            self.thisptr = new cudaMetropolis(
                rng,
                steps,
                scaling,
                acceptanceWeight,
                rejectionWeight,
                )
        
        def __dealloc__(self):
            # XXX: this will be freed by the base class
            #del self.thisptr
            pass

