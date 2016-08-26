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
wrapperSetupCurandState(curandState * const state, const int seed, const int Npoints)
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
wrapperMC(const float * const gpu_coor, const float * const gpu_radii, float * const gpu_points,
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
sascalc::cudaGV::
wrapperCuda_calcIq(const float * const gpu_coor, const float * const gpu_b, const int Natoms)
{
	// set the CUDA block dimenstions
    dim3 dim_grid(_NGV), dim_block(_Nq);
	// cudakernel launch configuration
	CALL_CUDA(cudaFuncSetCacheConfig(sascalc::cuda::cudaKernel_GV_calc, cudaFuncCachePreferL1));
	// launch the kernel
	sascalc::cuda::cudaKernel_GV_calc<<<dim_grid, dim_block>>>(gpu_coor, gpu_b, _gpu_gv, _gpu_q, _gpu_Is, Natoms, _Nq, _NGV);
	CALL_CUDA(cudaGetLastError());
}

/// @par Main functionality
/// wrap the cuda function by a C++ interface
/// @par CUDA threads layout
///- one thread corresponds to one sample
void
sascalc::cudaDebye::
wrapperCuda_calcIq(const float * const gpu_coor, const float * const gpu_b)
{
	// set the CUDA block dimenstions
    dim3 dim_grid((_Natoms-1)/BLOCKDIM+1), dim_block(BLOCKDIM);
	// cudakernel launch configuration
	CALL_CUDA(cudaFuncSetCacheConfig(sascalc::cuda::cudaKernel_Debye_calc, cudaFuncCachePreferShared));
	// launch the kernel
	sascalc::cuda::cudaKernel_Debye_calc<<<dim_grid, dim_block>>>(gpu_coor, gpu_b, _gpu_q, _gpu_Ia, _Natoms, _Nq);
	CALL_CUDA(cudaGetLastError());
}
