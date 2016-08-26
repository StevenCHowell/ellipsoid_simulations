

IF not IMPLEMENTATION:
    
    cdef extern from "altar/bayesian/SequentialAnnealing.h" namespace "altar::bayesian":
        
        cdef cppclass SequentialAnnealing(AnnealingMethod):
            
            SequentialAnnealing(size_t samples, size_t parameters, int rank) except +



cdef class XSequentialAnnealing(XAnnealingMethod):
    
    IF not IMPLEMENTATION:
        
        cdef SequentialAnnealing *thisptr   # hold a C++ instance which we're wrapping
    
    ELSE:
        
        def __cinit__(self, size_t samples, size_t parameters, int rank):
            self.thisptr = new SequentialAnnealing(samples, parameters, rank)
            self.annealingMethod = self.thisptr
        
        def __dealloc__(self):
            del self.thisptr
