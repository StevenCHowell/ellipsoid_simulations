//
// Hailiang Zhang
// california institute of technology
//

#ifndef H_CUDA_KERNELS_GINEMATIC
#define H_CUDA_KERNELS_GINEMATIC

// macros
#include <altar/utils/common.h>

#include <math.h>

#include <cuda.h>
#include <cuda_runtime.h>

#include <stdio.h>


// place everything in the local namespace
namespace altar {
    namespace bayesian {
        namespace cuda {
			/// \brief home of cuda kernel for @c cudaUpdateSample
			namespace cudaKernels_Metropolis {

				/// \brief Update @b M and <i>llk</i>'s for the accepted moves on GPU
				__global__ void cudaUpdateSample(const int Ns, const int NMparam, TYPE * const gR_select, const int * const rejects, int * gpnacc,
							TYPE * const gM, TYPE * const gllkprior, TYPE * const gllkdata, TYPE * const gllkpost,
							const TYPE * const gM_candidate, const TYPE * const gllkprior_candidate, const TYPE * const gllkdata_candidate, const TYPE * const gllkpost_candidate);

			} // of namespace cudaKernels_Metropolis
		} // of namespace cuda
    } // of namespace bayesian 
} // of namespace altar

#endif
