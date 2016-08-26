//
// Hailiang Zhang
// california institute of technology
//


#include "util_cuda.h"



//////////////////////////////////////////////////////////////
// get the flop of cublasDegemm
//////////////////////////////////////////////////////////////
long long
altar::bayesian::util::cuda::
get_flop_dgemm(long long m, long long n, long long k)
{
    return 2*2*m*n*k;
}



//////////////////////////////////////////////////////////////
// Set up the start cuda event time
//////////////////////////////////////////////////////////////
void
altar::bayesian::util::cuda::
init_cuda_event_time_diff(cudaEvent_t & start)
{
    CALL_CUDA( cudaEventRecord( start, 0 ) );
}


//////////////////////////////////////////////////////////////
// Set up the start cuda event time
//////////////////////////////////////////////////////////////
float
altar::bayesian::util::cuda::
get_cuda_event_time_diff(cudaEvent_t & start, cudaEvent_t & stop)
{
    float elapsedTime;
    CALL_CUDA( cudaEventRecord( stop, 0 ) );
    CALL_CUDA( cudaEventSynchronize( stop ) );
    CALL_CUDA( cudaEventElapsedTime( &elapsedTime,start, stop ) );
    return elapsedTime;
}
