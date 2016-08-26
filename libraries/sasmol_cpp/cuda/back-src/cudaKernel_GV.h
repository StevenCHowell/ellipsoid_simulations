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
        __global__ void cudaKernel_GV_calc(const float * const gpu_coor, const float * const gpu_b, const float * const gpu_gv, const float * const gpu_q, float * const gpu_I, const int Natoms, const int Nq, const int Ngv);
        __global__ void cudaKernel_GV_calc_subset(const float * const gpu_coor, const float * const gpu_b, const float * const gpu_gv, const float * const gpu_q, float * const gpu_I, const int Natoms, const int Natoms_subset, const int Nq, const int Ngv, const float scale);

        __global__ void cudaKernel_GV_calc_GVVV(const float * const gpu_coor, const float * const gpu_radius, const float * const gpu_gv, const float * const gpu_q, float * const gpu_I, const int Natoms, const int Nq, const int Ngv, const float sld_solvent);
        __global__ void cudaKernel_GV_calc_GVGaussian(const float * const gpu_coor, const float * const gpu_radius, const float * const gpu_gv, const float * const gpu_q, float * const gpu_I, const int Natoms, const float volume, const int Nq, const int Ngv, const float sld_solvent);

        __global__ void cudaKernel_GV_calc_Merge_ND(const float * const gpu_coor, const float * const gpu_b, const float * const gpu_radius, const float * const gpu_mc_coor, const float * const gpu_gv, const float * const gpu_q, float * const gpu_I, const int Natoms, const int N_mc_points, const float volume, const int Nq, const int Ngv, const float sld_solvent);
        __global__ void cudaKernel_GV_calc_Merge_NDsubset(const float * const gpu_coor, const float * const gpu_b, const float * const gpu_radius, const float * const gpu_mc_coor, const float * const gpu_gv, const float * const gpu_q, float * const gpu_I, const int Natoms, const int N_mc_points, const int N_mc_points_subset, const float volume, const int Nq, const int Ngv, const float sld_solvent);
        __global__ void cudaKernel_GV_calc_Merge_GVVV(const float * const gpu_coor, const float * const gpu_b, const float * const gpu_radius, const float * const gpu_gv, const float * const gpu_q, float * const gpu_I, const int Natoms, const int Nq, const int Ngv, const float sld_solvent, const float scale_volume);
        __global__ void cudaKernel_GV_calc_Merge_GVGaussian(const float * const gpu_coor, const float * const gpu_b, const float * const gpu_radius, const float * const gpu_gv, const float * const gpu_q, float * const gpu_I, const int Natoms, const float volume, const int Nq, const int Ngv, const float sld_solvent);
    }
}

#endif
