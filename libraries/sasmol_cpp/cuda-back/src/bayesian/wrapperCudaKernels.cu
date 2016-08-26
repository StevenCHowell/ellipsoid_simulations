//
// Hailiang Zhang
// california institute of technology
//

// My dependencies
#include <altar/utils/util.h>
#include <altar/utils/util_cuda.h>

#include <cuda.h>
#include <cuda_runtime.h>
#include <curand_kernel.h>
#include "cudaKernels.h"

#include <altar/bayesian/cudaMetropolis.h>

/// @par Main functionality
/// wrap the cudaUpdateSample function by a C++ interface <br>
/// @par CUDA threads layout
///- one thread corresponds to one sample
///- the number of threads per block is BLOCKDIM (defined in @c altar/utils/common.h)
/// @note see @c altar/bayesian/cudaKernels.cu for detailed parameter description
void
altar::bayesian::cuda::cudaMetropolis::
wrapperCudaUpdateSample(const int Ns, const int NMparam, const int * const & rejects, TYPE * const & gR_select, int * gpnacc,
						TYPE * const & gM, TYPE * const & gllkprior, TYPE * const & gllkdata, TYPE * const & gllkpost,
						const TYPE * const & gM_candidate, const TYPE * const & gllkprior_candidate, const TYPE * const & gllkdata_candidate, const TYPE * const & gllkpost_candidate) const
{
	// set the CUDA block dimenstions
    dim3 dim_grid((Ns-1)/BLOCKDIM+1), dim_block(BLOCKDIM);
	// cudakernel launch configuration
	CALL_CUDA(cudaFuncSetCacheConfig(altar::bayesian::cuda::cudaKernels_Metropolis::cudaUpdateSample, cudaFuncCachePreferL1)); 
	// launch the kernel
	altar::bayesian::cuda::cudaKernels_Metropolis::cudaUpdateSample<<<dim_grid, dim_block>>>(Ns, NMparam, gR_select, rejects, gpnacc,
					gM, gllkprior, gllkdata, gllkpost,
					gM_candidate, gllkprior_candidate, gllkdata_candidate, gllkpost_candidate);
	CALL_CUDA(cudaGetLastError());
}

