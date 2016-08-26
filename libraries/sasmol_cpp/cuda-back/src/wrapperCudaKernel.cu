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
wrapper_cudaMC_fill(const float * const gpu_coor, const float * const gpu_radii, float * const gpu_points,
                    const int Natoms, const int Npoints,
                    const float x0, const float y0, const float z0, const float xl, const float yl, const float zl,
                    curandState * const state)
{
	// set the CUDA block dimenstions
    dim3 dim_grid((Npoints-1)/BLOCKDIM+1), dim_block(BLOCKDIM);
	// cudakernel launch configuration
	CALL_CUDA(cudaFuncSetCacheConfig(cudaMC, cudaFuncCachePreferL1));
	// launch the kernel
	cudaMC<<<dim_grid, dim_block>>>(gpu_coor, gpu_radii, gpu_points, Natoms, Npoints, x0, y0, z0, xl, yl, zl, state);
	CALL_CUDA(cudaGetLastError());
}

/// @par Main functionality
/// wrap the cuda function by a C++ interface
/// @par CUDA threads layout
///- one thread corresponds to one sample
void
sascalc::cuda::
wrapper_cudaGV_calcIq(const float * const gpu_coor, const float * const gpu_b, const float * const gpu_gv, const float * const gpu_q, float * const gpu_Is, const int Natoms, const int Nq, const int NGV)
{
	// set the CUDA block dimenstions
    dim3 dim_grid(NGV), dim_block(Nq);
	// cudakernel launch configuration
	CALL_CUDA(cudaFuncSetCacheConfig(sascalc::cuda::cudaKernel_GV_calc, cudaFuncCachePreferL1));
	// launch the kernel
	sascalc::cuda::cudaKernel_GV_calc<<<dim_grid, dim_block>>>(gpu_coor, gpu_b, gpu_gv, gpu_q, gpu_Is, Natoms, Nq, NGV);
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
wrapper_cudaDebye_calcIq_CustomizedModel(const float * const gpu_coor, const float * const gpu_b, const float * const gpu_radius, const float * const gpu_q, float * const gpu_Ia, const int Natoms, const int Nq)
{
	// set the CUDA block dimenstions
    dim3 dim_grid((Natoms-1)/BLOCKDIM+1), dim_block(BLOCKDIM);
	// cudakernel launch configuration
	CALL_CUDA(cudaFuncSetCacheConfig(sascalc::cuda::cudaKernel_Debye_calc, cudaFuncCachePreferShared));
	// launch the kernel
	sascalc::cuda::cudaKernel_Debye_calc_CustomizedModel<<<dim_grid, dim_block>>>(gpu_coor, gpu_b, gpu_radius, gpu_q, gpu_Ia, Natoms, Nq);
	CALL_CUDA(cudaGetLastError());
}
