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
        __global__ void cudaKernel_GV_calc(const TYPE * const gpu_coor, const TYPE * const gpu_b, const TYPE * const gpu_gv, const TYPE * const gpu_q, TYPE * const gpu_I, const int Natoms, const int Nq, const int Ngv);
        __global__ void cudaKernel_GV_calc_subset(const TYPE * const gpu_coor, const TYPE * const gpu_b, const TYPE * const gpu_gv, const TYPE * const gpu_q, TYPE * const gpu_I, const int Natoms, const int Natoms_subset, const int Nq, const int Ngv, const TYPE scale);

        __global__ void cudaKernel_GV_calc_GVVV(const TYPE * const gpu_coor, const TYPE * const gpu_radius, const TYPE * const gpu_gv, const TYPE * const gpu_q, TYPE * const gpu_I, const int Natoms, const int Nq, const int Ngv, const TYPE sld_solvent);
        __global__ void cudaKernel_GV_calc_GVGaussian(const TYPE * const gpu_coor, const TYPE * const gpu_radius, const TYPE * const gpu_gv, const TYPE * const gpu_q, TYPE * const gpu_I, const int Natoms, const TYPE volume, const int Nq, const int Ngv, const TYPE sld_solvent);

        __global__ void cudaKernel_GV_calc_Merge_ND(const TYPE * const gpu_coor, const TYPE * const gpu_b, const TYPE * const gpu_radius, const TYPE * const gpu_mc_coor, const TYPE * const gpu_gv, const TYPE * const gpu_q, TYPE * const gpu_I, const int Natoms, const int N_mc_points, const TYPE volume, const int Nq, const int Ngv, const TYPE sld_solvent);
        __global__ void cudaKernel_GV_calc_Merge_NDsubset(const TYPE * const gpu_coor, const TYPE * const gpu_b, const TYPE * const gpu_radius, const TYPE * const gpu_mc_coor, const TYPE * const gpu_gv, const TYPE * const gpu_q, TYPE * const gpu_I, const int Natoms, const int N_mc_points, const int N_mc_points_subset, const TYPE volume, const int Nq, const int Ngv, const TYPE sld_solvent);
        __global__ void cudaKernel_GV_calc_Merge_GVVV(const TYPE * const gpu_coor, const TYPE * const gpu_b, const TYPE * const gpu_radius, const TYPE * const gpu_gv, const TYPE * const gpu_q, TYPE * const gpu_I, const int Natoms, const int Nq, const int Ngv, const TYPE sld_solvent, const TYPE scale_volume);
        __global__ void cudaKernel_GV_calc_Merge_GVGaussian(const TYPE * const gpu_coor, const TYPE * const gpu_b, const TYPE * const gpu_radius, const TYPE * const gpu_gv, const TYPE * const gpu_q, TYPE * const gpu_I, const int Natoms, const TYPE volume, const int Nq, const int Ngv, const TYPE sld_solvent);
    }
}

#endif
