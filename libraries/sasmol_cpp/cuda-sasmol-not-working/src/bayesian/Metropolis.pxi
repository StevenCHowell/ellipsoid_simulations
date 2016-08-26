

IF not IMPLEMENTATION:
    
    cdef extern from "altar/bayesian/Metropolis.h" namespace "altar::bayesian":
        
        cdef cppclass Metropolis:
           
            Metropolis() except +
             
            Metropolis(
                rng_t *,
                size_t steps,
                double scaling,
                double acceptanceWeight,
                double rejectionWeight,
                ) except +



cdef class XMetropolis:
    
    IF not IMPLEMENTATION:
        
        cdef Metropolis *thisptr   # hold a C++ instance which we're wrapping
        
        cdef object rngCapsule
    
    ELSE:
        
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
            self.thisptr = new Metropolis(
                rng,
                steps,
                scaling,
                acceptanceWeight,
                rejectionWeight,
                )
        
        def __dealloc__(self):
            del self.thisptr
