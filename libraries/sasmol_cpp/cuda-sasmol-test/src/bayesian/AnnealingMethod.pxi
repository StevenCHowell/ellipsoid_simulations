

IF not IMPLEMENTATION:
    
    cdef extern from "altar/bayesian/AnnealingMethod.h" namespace "altar::bayesian":
        
        cdef cppclass AnnealingMethod:
            
            pass



cdef class XAnnealingMethod:
    
    IF not IMPLEMENTATION:
        
        cdef AnnealingMethod *annealingMethod   # initialized by concrete subclasses
    
    ELSE:
        
        pass
