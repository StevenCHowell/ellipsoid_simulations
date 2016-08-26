//
// Hailiang Zhang
// NIST & UTK
//

#ifndef H_CUDA_KERNEL_GV
#define H_CUDA_KERNEL_GV

// macros
#include "cudaCommon.h"

#include <math.h>
#include <stdio.h>
#include <cuda.h>
#include <cuda_runtime.h>
#include <curand_kernel.h>


namespace sascalc
{
    namespace cuda
    {
        __global__ void cudaKernel_GV_calc(const float * const gpu_coor, const float * const gpu_b, const float * const gpu_gv, const float * const gpu_q, float * const gpu_I, const int Natoms, const int Nq, const int NGV);
    }
}

#endif
