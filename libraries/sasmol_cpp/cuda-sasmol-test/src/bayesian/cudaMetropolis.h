// -*- C++ -*-
//
// Hailiang Zhang
// california institute of technology
//

// code guard
#if !defined(altar_bayesian_cuda_cudaMetropolis_h)
#define altar_bayesian_cuda_cudaMetropolis_h

// macros
#include <altar/utils/common.h>

// place everything in the local namespace
namespace altar {
    namespace problem {
		namespace cuda {

			// forward declarations
			class cudaCoolingStep;
			class cudaProblem;
			
		} // of namespace cuda
    } // of namespace problem

    namespace bayesian {
		namespace cuda {

	        // forward declarations
        	class cudaMetropolis;

		} // of namespace cuda
    } // of namespace bayesian
} // of namespace altar

// externals
#include <iomanip>
#include <sys/stat.h>
#include <sstream>
#include <fstream>
#include <sys/mman.h>
#include <cmath>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#ifdef USE_DOUBLE 
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_statistics.h>
#else
#include <gsl/gsl_vector_float.h>
#include <gsl/gsl_matrix_float.h>
#include <gsl/gsl_statistics_float.h>
#endif
#include <cublas.h>
#include <cublas_v2.h>
#include <curand.h>

// my superclass 
#include <altar/bayesian/Metropolis.h>


// utilities
#include <altar/utils/util_cuda.h>

/// @brief Metropolis random walk on GPU
///
/// @par Primary contents
/// the number of chain steps (@i inherited)<br>
/// the scaling factor to be applied on @b Cs<br>
/// the empirical weights for acceptance/rejection<br>
/// the random number generator
///
/// @par Main functionalities
/// Random walk on GPU
///
/// @note
/// <i> Note to developers:</i><br>
/// cudaMetropolis has its own curand generator and cublas handle, primarily used for random sample generation on GPU
///

// declaration
class altar::bayesian::cuda::cudaMetropolis :
	public altar::bayesian::Metropolis
{
    // types
public:
    typedef altar::problem::cuda::cudaCoolingStep cuda_state_t;///< alias for cudaCoolingStep
    typedef altar::problem::cuda::cudaProblem cuda_problem_t;///< alias for cudaProblem

	// data
public:
	cublasHandle_t _cublas_handle;///< cublas handle
	curandGenerator_t _curand_gen;///< curand generator

    // interface
public:
	/// Random walk on GPU
    virtual stats_t sample(state_t &, cuda_state_t &, const cuda_problem_t &) const;

    // meta-methods
public:
	/// destructor
    virtual ~cudaMetropolis();
	/// constructor
    inline cudaMetropolis(
                      rng_t *, // the random number generator
                      size_t steps = 20, // the length of my Markov chains
                      TYPE scaling = .1, // factor by which the covariance matrix is scaled
                      TYPE acceptanceWeight = 8, // the relative weight of accepted samples
                      TYPE rejectionWeight = 1 // the relative weight of rejected samples
                      );

	// kernel wrappers                      
public:
	/// cuda kernel wrapper for @c cudaUpdateSample
	void wrapperCudaUpdateSample(const int Ns, const int NMparam, const int * const & rejects, TYPE * const & gR_select, int * gpnacc,
						TYPE * const & gM, TYPE * const & gllkprior, TYPE * const & gllkdata, TYPE * const & gllkpost,
						const TYPE * const & gM_candidate, const TYPE * const & gllkprior_candidate, const TYPE * const & gllkdata_candidate, const TYPE * const & gllkpost_candidate) const;

    // disallow
private:
	/// copy constructor disallowed
    inline cudaMetropolis(const cudaMetropolis &);
	/// assign constructor disallowed
    inline const cudaMetropolis & operator=(const cudaMetropolis &);
};

// get the inline definitions
#define altar_bayesian_cuda_cudaMetropolis_icc
#include "cudaMetropolis.icc"
#undef altar_bayesian_cuda_cudaMetropolis_icc

# endif
// end of file
