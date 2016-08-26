//
// Hailiang Zhang
// california institute of technology
//

#if !defined(altar_bayesian_utilcuda_h)
#define altar_bayesian_utilcuda_h

// My dependencies
#include <altar/utils/util_cuda.h>

#include <stdio.h>
#include <iostream>
#include <string.h>
#include <cuda.h>
#include <cuda_runtime.h>
#include <cuda_runtime_api.h>
#include <curand_kernel.h>
#include <cublas.h>
#include <cublas_v2.h>


// place everything in the local namespace
namespace altar {
    namespace bayesian {
		namespace util {
			namespace cuda {
				// get the flop of cublasDegemm
				long long get_flop_dgemm(long long m, long long n, long long k); 
				// Set up the start cuda event time
				void init_cuda_event_time_diff(cudaEvent_t & start);
				// Set up the start cuda event time
				float get_cuda_event_time_diff(cudaEvent_t & start, cudaEvent_t & stop);
				// print out the gpu data as a matrix
				// note that m is the leading index as Fortran style
				template <class DataType> void print_gpu_matrix(const DataType * const gpu, const int m, const int n, const int mp, const int np, const std::string message="", const std::string fmt="%16.3f")
				{
				    if (strlen(message.c_str())) std::cout<<message<<std::endl;
				    DataType *cpu=(DataType *)malloc(m*n*sizeof(DataType));
			    	CALL_CUDA(cudaMemcpy((void*)cpu, (void*)gpu, m*n*sizeof(DataType), cudaMemcpyDeviceToHost));
				    for (int i=0; i<np; i++)
				    {
				        for (int j=0; j<mp; j++) printf(fmt.c_str(),cpu[i*m+j]);
				        std::cout<<std::endl;
			    	}
				    std::cout<<std::endl;
				    free(cpu);
				}
			} // of namespace cuda
		} //namespace util
	} // namespace bayesian
} // namespace altar

#endif
