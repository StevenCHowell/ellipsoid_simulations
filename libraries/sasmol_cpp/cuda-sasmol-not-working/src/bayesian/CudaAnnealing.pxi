

IF not IMPLEMENTATION:
    
    cdef extern from "altar/bayesian/CudaAnnealing.h" namespace "altar::bayesian::cuda":
        
        cdef cppclass CudaAnnealing(AnnealingMethod):
            
            CudaAnnealing(size_t samples, size_t parameters, size_t ngpu, int rank) except +



cdef class XCudaAnnealing(XAnnealingMethod):
    
    IF not IMPLEMENTATION:
        
        cdef CudaAnnealing *thisptr   # hold a C++ instance which we're wrapping
    
    ELSE:
        
        def __cinit__(self, size_t samples, size_t parameters, size_t ngpu, int rank):
            self.thisptr = new CudaAnnealing(samples, parameters, ngpu, rank)
            self.annealingMethod = self.thisptr
        
        def __dealloc__(self):
            del self.thisptr
