// -*- C++ -*-
//
// Hailiang Zhang
// california institute of technology
//

// code guard
#if !defined(altar_bayesian_cuda_CudaAnnealing_h)
#define altar_bayesian_cuda_CudaAnnealing_h

// macros
#include <altar/utils/common.h>

// place everything in the local namespace
namespace altar {
    namespace problem {
        class cudaCoolingStep;
    }
    namespace bayesian {
        namespace cuda {

            // forward declarations
            class CudaAnnealing;

        } // of namespace cuda
    } // of namespace bayesian
} // of namespace altar

// my superclass
#include <altar/bayesian/AnnealingMethod.h>

// my dependencies
#include "cudaMetropolis.h"
#include <altar/bayesian/Annealer.h>
#include <altar/problem/CoolingStep.h>
#include <altar/problem/cudaCoolingStep.h>
#include <altar/utils/util.h>
#include <altar/utils/util_cuda.h>

// external
#include <string>

/// @brief Annealing events on GPU
///
/// @par Primary contents
/// @c cuda_state: the place where @b M, <i>LLK</i>'s and @b Cs resides on GPU
///
/// @par Main functionalities
/// start the annealing process from scratch<br>
/// restart the annealing process from a previous run<br>
/// resample<br>
/// status report
///
/// @note
/// <i> Note to developers:</i><br>
/// none
/// @todo
/// The MPI stuff needs to be removed
///

// base class for the various annealing methods
class altar::bayesian::cuda::CudaAnnealing :
    public altar::bayesian::AnnealingMethod
{

    // types
public:
    typedef cudaMetropolis cuda_sampler_t;///< alias for cudaMetropolis
    typedef altar::problem::cuda::cudaCoolingStep cuda_state_t;///< alias for cudaCoolingStep
    typedef altar::problem::cuda::cudaProblem cuda_problem_t;///< alias for cudaProblem

    // data
protected:
    cuda_state_t cuda_state;///<the place where @b M, <i>LLK</i>'s and @b Cs resides on GPU

    // interface
public:
    /// start the annealing process from scratch
    virtual AnnealingMethod & start(annealer_t &, const std::string &);
    /// restart the annealing process from a previous run
    virtual AnnealingMethod & restart(annealer_t &, const size_t iteration);
    /// the annealing steps
    virtual stats_t resample(annealer_t &);

    // meta-methods
public:
    /// destructor
    virtual ~CudaAnnealing();
    /// constructor
    inline CudaAnnealing(const size_t samples, const size_t parameters, const size_t ngpu, const unsigned short rank=0);

    // implementation details
protected:
    /// status report engine
    virtual AnnealingMethod & _report(annealer_t &);

    // disallow
private:
    /// copy constructor disallowed
    inline CudaAnnealing(const CudaAnnealing &);
    /// assign constructor disallowed
    inline const CudaAnnealing & operator=(const CudaAnnealing &);
};

// get the inline definitions
#define altar_bayesian_cuda_CudaAnnealing_icc
#include "CudaAnnealing.icc"
#undef altar_bayesian_cuda_CudaAnnealing_icc

# endif
// end of file
