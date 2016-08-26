

IF not IMPLEMENTATION:
    
    cdef extern from "altar/bayesian/MPIAnnealing.h" namespace "altar::bayesian":
        
        cdef cppclass MPIAnnealing(AnnealingMethod):
            
            MPIAnnealing(
                communicator_t communicator,
                size_t samples,
                size_t parameters,
                AnnealingMethod worker,
                ) except +



cdef class XMPIAnnealing(XAnnealingMethod):
    
    IF not IMPLEMENTATION:
        
        cdef MPIAnnealing *thisptr   # hold a C++ instance which we're wrapping
        
        cdef object communicatorCapsule
    
    ELSE:
        
        def __cinit__(self,
                      communicatorCapsule,
                      size_t samples,
                      size_t parameters,
                      XAnnealingMethod worker,
                      ):
            # retain my children
            self.communicatorCapsule = communicatorCapsule
            cdef communicator_t *communicator
            communicator = <communicator_t *>PyCapsule_GetPointer(self.communicatorCapsule, communicator_capsule_name)
            self.thisptr = new MPIAnnealing(
                communicator[0],
                samples,
                parameters,
                worker.annealingMethod[0],
                )
            self.annealingMethod = self.thisptr
        
        def __dealloc__(self):
            del self.thisptr
