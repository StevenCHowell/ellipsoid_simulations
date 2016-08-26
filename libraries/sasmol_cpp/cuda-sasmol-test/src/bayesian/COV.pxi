

IF not IMPLEMENTATION:
    
    cdef extern from "altar/bayesian/COV.h" namespace "altar::bayesian":
        
        cdef cppclass COV:
            
            COV(rng_t *rng,
                double tolerance,
                size_t maxIterations,
                double target,
                ) except +



cdef class XCOV:
    
    IF not IMPLEMENTATION:
        
        cdef COV *thisptr   # hold a C++ instance which we're wrapping
        
        cdef object rngCapsule
    
    ELSE:
        
        def __cinit__(self,
                      rngCapsule,
                      double tolerance,
                      size_t maxIterations,
                      double target,
                      ):
            # retain my children
            self.rngCapsule = rngCapsule
            cdef rng_t *rng
            rng = <rng_t *>PyCapsule_GetPointer(self.rngCapsule, rng_capsule_name)
            self.thisptr = new COV(rng, tolerance, maxIterations, target)
        
        def __dealloc__(self):
            del self.thisptr
