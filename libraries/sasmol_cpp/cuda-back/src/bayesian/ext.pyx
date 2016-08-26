
DEF IMPLEMENTATION = True

include "ext.pxi"

IF WITH_CUDA:
    # XXX: ugly hack which doesn't even belong here
    cdef extern from "altar/utils/util_cuda.h":
        int cudaSetDevice(int n)
        void CALL_CUDA(int err)
    def setCudaDevice(n):
        CALL_CUDA(cudaSetDevice(n))

