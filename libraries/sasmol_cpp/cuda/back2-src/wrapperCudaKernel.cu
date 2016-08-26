//
// Hailiang Zhang
// NIST & UTK
//

#include "wrapperCudaKernel.h"

#include "cudaKernel_MC.h"
#include "cudaKernel_GV.h"
#include "cudaKernel_Debye.h"

#include "cudaUtil.h"

#include <iostream>

#include <cuda.h>
#include <cuda_runtime.h>
#include <curand_kernel.h>


/// @par Main functionality
/// wrap the cuda function by a C++ interface
/// @par CUDA threads layout
///- one thread corresponds to one sample
void
sascalc::cuda::
wrapper_cudaMC_setupCurandState(curandState * const state, const int seed, const int Npoints)
{
	// set the CUDA block dimenstions
    dim3 dim_grid((Npoints-1)/BLOCKDIM+1), dim_block(BLOCKDIM);
	// cudakernel launch configuration
	CALL_CUDA(cudaFuncSetCacheConfig(setupCurandKernel, cudaFuncCachePreferL1));
	// launch the kernel
	setupCurandKernel<<<dim_grid, dim_block>>>(state, seed, Npoints);
	CALL_CUDA(cudaGetLastError());
}

/// @par Main functionality
/// wrap the cuda function by a C++ interface
/// @par CUDA threads layout
///- one thread corresponds to one sample
void
sascalc::cuda::
wrapper_cudaMC_setupGrid(const float * const gpu_coor, const float * const gpu_radii, int * const gpu_grid, const int Natoms,
                         const float x0, const float y0, const float z0, const float delta, const int Nx, const int Ny, const int Nz)
{
	// set the CUDA block dimenstions
    dim3 dim_grid((Natoms-1)/BLOCKDIM+1), dim_block(BLOCKDIM);
	// cudakernel launch configuration
	CALL_CUDA(cudaFuncSetCacheConfig(cudaSetupGrid, cudaFuncCachePreferL1));
	// launch the kernel
	cudaSetupGrid<<<dim_grid, dim_block>>>(gpu_coor, gpu_radii, gpu_grid, Natoms, x0, y0, z0, delta, Nx, Ny, Nz);
	CALL_CUDA(cudaGetLastError());
}

/// @par Main functionality
/// wrap the cuda function by a C++ interface
/// @par CUDA threads layout
///- one thread corresponds to one sample
void
sascalc::cuda::
wrapper_cudaMC_fill(const float * const gpu_coor, const float * const gpu_radii, const int * const gpu_grid, float * const gpu_points,
                    const int Natoms, const int Npoints, int * const gpu_Npoints_tried,
                    const float x0, const float y0, const float z0, const float delta, const int Nx, const int Ny, const int Nz,
                    curandState * const state)
{
	// set the CUDA block dimenstions
    dim3 dim_grid((Npoints-1)/BLOCKDIM+1), dim_block(BLOCKDIM);
	// cudakernel launch configuration
	CALL_CUDA(cudaFuncSetCacheConfig(cudaMC, cudaFuncCachePreferL1));
	// launch the kernel
	cudaMC<<<dim_grid, dim_block>>>(gpu_coor, gpu_radii, gpu_grid, gpu_points, Natoms, Npoints, gpu_Npoints_tried, x0, y0, z0, delta, Nx, Ny, Nz, state);
	CALL_CUDA(cudaGetLastError());
}

/// @par Main functionality
/// wrap the cuda function by a C++ interface
/// @par CUDA threads layout
///- one thread corresponds to one sample
void
sascalc::cuda::
wrapper_cudaGV_calcIq(const float * const gpu_coor, const float * const gpu_b, const float * const gpu_gv, const float * const gpu_q, float * const gpu_Is, const int Natoms, const int Nq, const int Ngv)
{
	// set the CUDA block dimenstions
    dim3 dim_grid(Ngv), dim_block(Nq);
	// cudakernel launch configuration
	CALL_CUDA(cudaFuncSetCacheConfig(sascalc::cuda::cudaKernel_GV_calc, cudaFuncCachePreferL1));
	// launch the kernel
	sascalc::cuda::cudaKernel_GV_calc<<<dim_grid, dim_block>>>(gpu_coor, gpu_b, gpu_gv, gpu_q, gpu_Is, Natoms, Nq, Ngv);
	CALL_CUDA(cudaGetLastError());
}

/// @par Main functionality
/// wrap the cuda function by a C++ interface
/// @par CUDA threads layout
///- one thread corresponds to one sample
void
sascalc::cuda::
wrapper_cudaGV_calcIq_subset(const float * const gpu_coor, const float * const gpu_b, const float * const gpu_gv, const float * const gpu_q, float * const gpu_Is, const int Natoms, const int Natoms_subset, const int Nq, const int Ngv, const float scale)
{
	// set the CUDA block dimenstions
    dim3 dim_grid(Ngv), dim_block(Nq);
	// cudakernel launch configuration
	CALL_CUDA(cudaFuncSetCacheConfig(sascalc::cuda::cudaKernel_GV_calc_subset, cudaFuncCachePreferL1));
	// launch the kernel
	sascalc::cuda::cudaKernel_GV_calc_subset<<<dim_grid, dim_block>>>(gpu_coor, gpu_b, gpu_gv, gpu_q, gpu_Is, Natoms, Natoms_subset, Nq, Ngv, scale);
	CALL_CUDA(cudaGetLastError());
}

/// @par Main functionality
/// wrap the cuda function by a C++ interface
/// @par CUDA threads layout
///- one thread corresponds to one sample
void
sascalc::cuda::
wrapper_cudaGV_calcIq_GVVV(const float * const gpu_coor, const float * const gpu_radius, const float * const gpu_gv, const float * const gpu_q, float * const gpu_Is, const int Natoms, const int Nq, const int Ngv, const float sld_solvent)
{
	// set the CUDA block dimenstions
    dim3 dim_grid(Ngv), dim_block(Nq);
	// cudakernel launch configuration
	CALL_CUDA(cudaFuncSetCacheConfig(sascalc::cuda::cudaKernel_GV_calc_GVVV, cudaFuncCachePreferL1));
	// launch the kernel
	sascalc::cuda::cudaKernel_GV_calc_GVVV<<<dim_grid, dim_block>>>(gpu_coor, gpu_radius, gpu_gv, gpu_q, gpu_Is, Natoms, Nq, Ngv, sld_solvent);
	CALL_CUDA(cudaGetLastError());
}

/// @par Main functionality
/// wrap the cuda function by a C++ interface
/// @par CUDA threads layout
///- one thread corresponds to one sample
void
sascalc::cuda::
wrapper_cudaGV_calcIq_GVGaussian(const float * const gpu_coor, const float * const gpu_radius, const float * const gpu_gv, const float * const gpu_q, float * const gpu_Is, const int Natoms, const float volume, const int Nq, const int Ngv, const float sld_solvent)
{
	// set the CUDA block dimenstions
    dim3 dim_grid(Ngv), dim_block(Nq);
	// cudakernel launch configuration
	CALL_CUDA(cudaFuncSetCacheConfig(sascalc::cuda::cudaKernel_GV_calc_GVGaussian, cudaFuncCachePreferL1));
	// launch the kernel
	sascalc::cuda::cudaKernel_GV_calc_GVGaussian<<<dim_grid, dim_block>>>(gpu_coor, gpu_radius, gpu_gv, gpu_q, gpu_Is, Natoms, volume, Nq, Ngv, sld_solvent);
	CALL_CUDA(cudaGetLastError());
}

/// @par Main functionality
/// wrap the cuda function by a C++ interface
/// @par CUDA threads layout
///- one thread corresponds to one sample
void
sascalc::cuda::
wrapper_cudaGV_calcIq_Merge_ND(const float * const gpu_coor, const float * const gpu_b, const float * const gpu_radius, const float * const gpu_mc_coor, const float * const gpu_gv, const float * const gpu_q, float * const gpu_Is, const int Natoms, const int N_mc_points, const float volume, const int Nq, const int Ngv, const float sld_solvent)
{
	// set the CUDA block dimenstions
    dim3 dim_grid(Ngv), dim_block(Nq);
	// cudakernel launch configuration
	CALL_CUDA(cudaFuncSetCacheConfig(sascalc::cuda::cudaKernel_GV_calc_Merge_ND, cudaFuncCachePreferL1));
	// launch the kernel
	sascalc::cuda::cudaKernel_GV_calc_Merge_ND<<<dim_grid, dim_block>>>(gpu_coor, gpu_b, gpu_radius, gpu_mc_coor, gpu_gv, gpu_q, gpu_Is, Natoms, N_mc_points, volume, Nq, Ngv, sld_solvent);
	CALL_CUDA(cudaGetLastError());
}

/// @par Main functionality
/// wrap the cuda function by a C++ interface
/// @par CUDA threads layout
///- one thread corresponds to one sample
void
sascalc::cuda::
wrapper_cudaGV_calcIq_Merge_GVVV(const float * const gpu_coor, const float * const gpu_b, const float * const gpu_radius, const float * const gpu_gv, const float * const gpu_q, float * const gpu_Is, const int Natoms, const int Nq, const int Ngv, const float sld_solvent)
{
	// set the CUDA block dimenstions
    dim3 dim_grid(Ngv), dim_block(Nq);
	// cudakernel launch configuration
	CALL_CUDA(cudaFuncSetCacheConfig(sascalc::cuda::cudaKernel_GV_calc_Merge_GVVV, cudaFuncCachePreferL1));
	// launch the kernel
	sascalc::cuda::cudaKernel_GV_calc_Merge_GVVV<<<dim_grid, dim_block>>>(gpu_coor, gpu_b, gpu_radius, gpu_gv, gpu_q, gpu_Is, Natoms, Nq, Ngv, sld_solvent);
	CALL_CUDA(cudaGetLastError());
}

/// @par Main functionality
/// wrap the cuda function by a C++ interface
/// @par CUDA threads layout
///- one thread corresponds to one sample
void
sascalc::cuda::
wrapper_cudaGV_calcIq_Merge_GVGaussian(const float * const gpu_coor, const float * const gpu_b, const float * const gpu_radius, const float * const gpu_gv, const float * const gpu_q, float * const gpu_Is, const int Natoms, const float volume, const int Nq, const int Ngv, const float sld_solvent)
{
	// set the CUDA block dimenstions
    dim3 dim_grid(Ngv), dim_block(Nq);
	// cudakernel launch configuration
	CALL_CUDA(cudaFuncSetCacheConfig(sascalc::cuda::cudaKernel_GV_calc_Merge_GVGaussian, cudaFuncCachePreferL1));
	// launch the kernel
	sascalc::cuda::cudaKernel_GV_calc_Merge_GVGaussian<<<dim_grid, dim_block>>>(gpu_coor, gpu_b, gpu_radius, gpu_gv, gpu_q, gpu_Is, Natoms, volume, Nq, Ngv, sld_solvent);
	CALL_CUDA(cudaGetLastError());
}

/// @par Main functionality
/// wrap the cuda function by a C++ interface
/// @par CUDA threads layout
///- one thread corresponds to one sample
void
sascalc::cuda::
wrapper_cudaGV_calcIq_Merge_NDsubset(const float * const gpu_coor, const float * const gpu_b, const float * const gpu_radius, const float * const gpu_mc_coor, const float * const gpu_gv, const float * const gpu_q, float * const gpu_Is, const int Natoms, const int N_mc_points, const int N_mc_points_subset, const float volume, const int Nq, const int Ngv, const float sld_solvent)
{
	// set the CUDA block dimenstions
    dim3 dim_grid(Ngv), dim_block(Nq);
	// cudakernel launch configuration
	CALL_CUDA(cudaFuncSetCacheConfig(sascalc::cuda::cudaKernel_GV_calc_Merge_NDsubset, cudaFuncCachePreferL1));
	// launch the kernel
	sascalc::cuda::cudaKernel_GV_calc_Merge_NDsubset<<<dim_grid, dim_block>>>(gpu_coor, gpu_b, gpu_radius, gpu_mc_coor, gpu_gv, gpu_q, gpu_Is, Natoms, N_mc_points, N_mc_points_subset, volume, Nq, Ngv, sld_solvent);
	CALL_CUDA(cudaGetLastError());
}

/// @par Main functionality
/// wrap the cuda function by a C++ interface
/// @par CUDA threads layout
///- one thread corresponds to one sample
void
sascalc::cuda::
wrapper_cudaDebye_calcIq(const float * const gpu_coor, const float * const gpu_b, const float * const gpu_q, float * const gpu_Ia, const int Natoms, const int Nq)
{
	// set the CUDA block dimenstions
    dim3 dim_grid((Natoms-1)/BLOCKDIM+1), dim_block(BLOCKDIM);
	// cudakernel launch configuration
	CALL_CUDA(cudaFuncSetCacheConfig(sascalc::cuda::cudaKernel_Debye_calc, cudaFuncCachePreferShared));
	// launch the kernel
	sascalc::cuda::cudaKernel_Debye_calc<<<dim_grid, dim_block>>>(gpu_coor, gpu_b, gpu_q, gpu_Ia, Natoms, Nq);
	CALL_CUDA(cudaGetLastError());
}

/// @par Main functionality
/// wrap the cuda function by a C++ interface
/// @par CUDA threads layout
///- one thread corresponds to one sample
void
sascalc::cuda::
wrapper_cudaDebye_calcIq_subset(const float * const gpu_coor, const float * const gpu_b, const float * const gpu_q, float * const gpu_Ia, const int Natoms, const int Natoms_subset, const int Nq, const float scale)
{
	// set the CUDA block dimenstions
    dim3 dim_grid((Natoms_subset-1)/BLOCKDIM+1), dim_block(BLOCKDIM);
	// cudakernel launch configuration
	CALL_CUDA(cudaFuncSetCacheConfig(sascalc::cuda::cudaKernel_Debye_calc_subset, cudaFuncCachePreferShared));
	// launch the kernel
	sascalc::cuda::cudaKernel_Debye_calc_subset<<<dim_grid, dim_block>>>(gpu_coor, gpu_b, gpu_q, gpu_Ia, Natoms, Natoms_subset, Nq, scale);
	CALL_CUDA(cudaGetLastError());
}

/// @par Main functionality
/// wrap the cuda function by a C++ interface
/// @par CUDA threads layout
///- one thread corresponds to one sample
void
sascalc::cuda::
wrapper_cudaDebye_calcIq_CustomizedModel(const float * const gpu_coor, const float * const gpu_b, const float * const gpu_radius, const float * const gpu_q, float * const gpu_Ia, const int Natoms, const int Nq)
{
	// set the CUDA block dimenstions
    dim3 dim_grid((Natoms-1)/BLOCKDIM+1), dim_block(BLOCKDIM);
	// cudakernel launch configuration
	CALL_CUDA(cudaFuncSetCacheConfig(sascalc::cuda::cudaKernel_Debye_calc_CustomizedModel, cudaFuncCachePreferShared));
	// launch the kernel
	sascalc::cuda::cudaKernel_Debye_calc_CustomizedModel<<<dim_grid, dim_block>>>(gpu_coor, gpu_b, gpu_radius, gpu_q, gpu_Ia, Natoms, Nq);
	CALL_CUDA(cudaGetLastError());
}
