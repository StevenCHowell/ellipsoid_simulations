//
// Hailiang Zhang
// NIST & UTK
//

#ifndef H_CUDA_KERNEL_DEBYE
#define H_CUDA_KERNEL_DEBYE

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
        __global__ void cudaKernel_Debye_calc(const float * const gpu_coor, const float * const gpu_b, const float * const gpu_q, float * const gpu_Ia, const int Natoms, const int Nq);
        __global__ void cudaKernel_Debye_calc_CustomizedModel(const float * const gpu_coor, const float * const gpu_b, const float * const gpu_radius, const float * const gpu_q, float * const gpu_Ia, const int Natoms, const int Nq);
    }
}


#endif
