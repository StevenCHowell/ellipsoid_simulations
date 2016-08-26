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
        __global__ void cudaKernel_Debye_calc(const TYPE * const gpu_coor, const TYPE * const gpu_b, const TYPE * const gpu_q, TYPE * const gpu_Ia, const int Natoms, const int Nq);
        __global__ void cudaKernel_Debye_calc_subset(const TYPE * const gpu_coor, const TYPE * const gpu_b, const TYPE * const gpu_q, TYPE * const gpu_Ia, const int Natoms, const int Natoms_subset, const int Nq, const TYPE scale);
        __global__ void cudaKernel_Debye_calc_CustomizedModel(const TYPE * const gpu_coor, const TYPE * const gpu_b, const TYPE * const gpu_radius, const TYPE * const gpu_q, TYPE * const gpu_Ia, const int Natoms, const int Nq);
    }
}


#endif
