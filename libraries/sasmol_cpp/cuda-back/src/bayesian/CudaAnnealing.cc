// -*- C++ -*-
// 
// Hailiang Zhang
// california institute of technology
// 


// for the build system
#include <portinfo>

// get my declarations
#include "CudaAnnealing.h"
#include <altar/problem/cudaProblem.h>
#include <altar/problem/cudaCoolingStep.h>

#include "rank.h"

/// @par Main functionality
/// Start the annealing process from scratch<br>
///- sample is initialized<br>
///- likelihoods are calculated
/// @param [in, out] annealer Host of the problem entity where priors, models and data reside
/// @param [in] fprior File name of the to-be-read prior
altar::bayesian::AnnealingMethod &
altar::bayesian::cuda::CudaAnnealing::
start(annealer_t & annealer, const std::string & fprior)
{
    // get the problem
    cuda_problem_t & cuda_problem = static_cast<cuda_problem_t &>(annealer.problem());

    // initialize the prior/model-ensemble (prior/data-LLK and rand-seeds
    cuda_problem.initialEnsemble(cuda_state.samples(), altar::utils::random_seed());

    // get the worker rank
    int rank = my_rank();

    // ask it to initialize/read my sample set
    if (!fprior.length()) cuda_problem.initialSample(cuda_state, rank);
    else cuda_problem.readSample(cuda_state, fprior, rank);

    // and compute the log likelihoods
    cuda_problem.likelihoods(cuda_state);

    // and send GPU data to CPU
    cuda_state.setLLKs();
    cuda_state.collectGPU(_state);

    // print-outs
    if (rank==0)
    {
        cuda_problem.show(cuda_state.samples());
        cuda_state.show();
    }

    // all done
    return AnnealingMethod::start(annealer, fprior);
}


/// @par Main functionality
/// Start the annealing process from a previous run<br>
/// @param [in, out] annealer Host of the problem entity where priors, models and data reside
/// @param [in] iteration The iteration number to be restarted
altar::bayesian::AnnealingMethod &
altar::bayesian::cuda::CudaAnnealing::
restart(annealer_t & annealer, const size_t iteration)
{
    // get the problem
    cuda_problem_t & cuda_problem = static_cast<cuda_problem_t &>(annealer.problem());

    // initialize the prior/model-ensemble (prior/data-LLK and rand-seeds
    cuda_problem.initialEnsemble(cuda_state.samples(), altar::utils::random_seed());

    // ask it to read the restart samples, llks and others
    cuda_state.initialRestart(iteration);

    // and send GPU data to CPU
    //cuda_state.setLLKs(); // do not need this because LLKs is already in regular place instead of candidate
    cuda_state.collectGPU(_state);

    // print-outs
    int rank = my_rank();
    if (rank==0)
    {
        cuda_problem.show(cuda_state.samples());
        cuda_state.show();
    }

    // all done
    return AnnealingMethod::restart(annealer, iteration);
}

/// @par Main functionality
/// Re-sample the posterior distribution
/// @param [in, out] annealer Host of the problem entity where priors, models and data reside
altar::bayesian::AnnealingMethod::stats_t
altar::bayesian::cuda::CudaAnnealing::
resample(annealer_t & annealer)
{
    // get the problem
    cuda_problem_t & cuda_problem = static_cast<cuda_problem_t &>(annealer.problem());
    // and the sampler
    cuda_sampler_t & cuda_sampler = static_cast<cuda_sampler_t &>(annealer.sampler());

    // ask it to sampler the posterior pdf and return the sampling statistics
    return cuda_sampler.sample(_state, cuda_state, cuda_problem);
}


/// @par Main functionality
/// Destructor
altar::bayesian::cuda::CudaAnnealing::
~CudaAnnealing() 
{}


/// @par Main functionality
/// Report progress
/// @param [in] annealer Host of the problem entity where priors, models and data reside
altar::bayesian::AnnealingMethod &
altar::bayesian::cuda::CudaAnnealing::
_report(annealer_t & annealer)
{
    // chain up
    AnnealingMethod::_report(annealer);

    // show me
    _info
        // the current temperature
        << "           beta: " << _state.beta() << pyre::journal::newline;

    const char * indent = "                 ";
    
    // show me {theta}
    state_t::matrix_t * m = _state.theta();
    _info
        << "          theta: (" << m->size1 
        << " samples) x (" << m->size2 << " parameters)" << pyre::journal::newline;
    // for small enough values, show me contents
    if (m->size1 < 10 && m->size2 < 10) {
        // for each row
        for (size_t row=0; row<m->size1; ++row) {
            _info << indent << (row == 0 ? "[[" : "  ");
            for (size_t column=0; column<m->size2; ++column) {
                _info << ' ' << gsl_matrix_TYPE_get(m, row, column);
            } 
            _info << (row == m->size1-1 ? "]]" : "") << pyre::journal::newline;
        }
    }

    // show me the prior
    state_t::vector_t * v = _state.prior();
    _info
        << "          prior: " << v->size << " samples" << pyre::journal::newline;
    // for small enough values, show me contents
    if (v->size < 10) {
        _info << indent << "[";
        for (size_t index=0; index<v->size; ++index) {
            _info << ' ' << gsl_vector_TYPE_get(v, index);
        } 
        _info << " ]" << pyre::journal::newline;
    }


    // show me the data
    v = _state.data();
    _info
        << "           data: " << v->size << " samples" << pyre::journal::newline;
    // for small enough values, show me contents
    if (v->size < 10) {
        _info << indent << "[";
        for (size_t index=0; index<v->size; ++index) {
            _info << ' ' << gsl_vector_TYPE_get(v, index);
        } 
        _info << " ]" << pyre::journal::newline;
    }

    // show me the posterior
    v = _state.posterior();
    _info
        << "      posterior: " << v->size << " samples" << pyre::journal::newline;
    // for small enough values, show me contents
    if (v->size < 10) {
        _info << indent << "[";
        for (size_t index=0; index<v->size; ++index) {
            _info << ' ' << gsl_vector_TYPE_get(v, index);
        } 
        _info << " ]" << pyre::journal::newline;
    }

    // show me {sigma}
    m = _state.sigma();
    _info
        << "          sigma: (" << m->size1 
        << " parameters) x (" << m->size2 << " parameters)" << pyre::journal::newline;
    // for small enough values, show me contents
    if (m->size1 < 10 && m->size2 < 10) {
        // for each row
        for (size_t row=0; row<m->size1; ++row) {
            _info << indent << (row == 0 ? "[[" : "  ");
            for (size_t column=0; column<m->size2; ++column) {
                _info << ' ' << gsl_matrix_TYPE_get(m, row, column);
            } 
            _info << (row == m->size1-1 ? "]]" : "") << pyre::journal::newline;
        }
    }

    // all done
    return *this;
}


// end of file
