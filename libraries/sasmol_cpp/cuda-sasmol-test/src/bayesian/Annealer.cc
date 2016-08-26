// -*- C++ -*-
// 
// michael a.g. aïvázis
// Hailiang Zhang
// california institute of technology
// (c) 2010-2013 all rights reserved
// 


// for the build system
#include <portinfo>
#include <pyre/journal.h>

// get my declarations
#include "Annealer.h"

// my dependencies
#include "COV.h"
#include "Metropolis.h"
#include "AnnealingMethod.h"
#include "SequentialAnnealing.h"

#include <altar/problem/CoolingStep.h>

#include "rank.h"

/// @par Main functionality
/// Sample the posterior probability density function of the given problem
/// @param [in] restart The restart iteration number
/// @param [in] wfreq The frequency to write the samples to the external file
/// @param [in] fprior File name of the to-be-read prior
int
altar::bayesian::Annealer::
posterior(const size_t restart, const size_t wfreq, const std::string & fprior)
{
    // make it all look prettier
    method_t & worker = _worker;
    // start from zero
    if (!restart) worker.start(*this, fprior);
    // restart from a previous run
    else worker.restart(*this, restart);

    // Write the initial samples and statistics
    stats_t stats1 = stats_t(0.0, 0.0);
    worker.write_out(*this, stats1, true, restart); 
    
    size_t niter = 0;
    // iterate until beta is sufficiently close to 1
    stats_t stats;

    // profiling preparation
    struct timeval t1, t2;
    bool P_altar_profiling = false;
    if (getenv("ALTAR_PROFILING") && strcmp(getenv("ALTAR_PROFILING"),"true")==0) P_altar_profiling = true;
    if (P_altar_profiling) gettimeofday(&t1,NULL);
    // get the rank
    int rank = my_rank();

    // remove all the exiiting profiling data
    std::string cuda_statistics_dir = "cuda_statistics/"; // it's important that directory name is postfixed with a slash
    if ( (rank==0) && opendir(cuda_statistics_dir.c_str()) )
    {
        system( (std::string("rm -f ")+cuda_statistics_dir+std::string("profile_rank_*")).c_str() );
    }

    // iterate
    while (worker.beta() + _tolerance < 1) {
        // Update Model Data (e.g., Cp stuff)
        worker.updateModelData(*this);
        // show me progress
        worker.report(*this);
        if (P_altar_profiling) altar::utils::print_profiling(rank,t1, t2, "report");
        // compute an acceptable temperature increment
        worker.cool(*this);
        if (P_altar_profiling) altar::utils::print_profiling(rank,t1, t2, "cooling");
        // re-sample
        stats = worker.resample(*this);
        if (P_altar_profiling) altar::utils::print_profiling(rank,t1, t2, "resample");
        // equilibrate
        worker.equilibrate(*this, stats);
        if (P_altar_profiling) altar::utils::print_profiling(rank,t1, t2, "equilibrate");
        // write out progress
        ++niter;
        worker.write_out(*this, stats, (wfreq>0 && niter%wfreq==0));
        if (P_altar_profiling) altar::utils::print_profiling(rank,t1, t2, "write_out");
    }

    if (worker.updateModelData(*this))
    {
        //printf("FINAL UPDATE\n");
        worker.report(*this);
        worker.iterate();
        worker.write_out(*this, stats, (wfreq>0 && niter%wfreq==0));
        stats = worker.resample(*this);
        //worker.equilibrate(*this, stats);
        ++niter;
    }

    // show me progress
    worker.report(*this);
    // write the final samples and statistics
    worker.write_out(*this, stats, true);
    // all done
    worker.finish(*this);
    
    // indicate success
    return 0;
}

/// @par Main functionality
/// Report progress
altar::bayesian::Annealer &
altar::bayesian::Annealer::
report()
{
    // ask my worker
    _worker.report(*this);
    // all done
    return *this;
}

/// @par Main functionality
/// Get the number of parameters
/// @return The number of parameters
size_t
altar::bayesian::Annealer::
parameters() const
{
    return _worker.state().parameters();
}


/// @par Main functionality
/// destructor
altar::bayesian::Annealer::
~Annealer() 
{}

/// @par Main functionality
/// constructor
/// @param [in, out] problem The problem entity (including priors, models, data)
/// @param [in, out] scheduler Class responsible for cooling between beta steps (alias for <i>COV</i>)
/// @param [in] sampler Class responsible for random walk during a beta step (alias for <i>Metropolis</i>)
/// @param [in] worker Worker to implement the annealing (alias for <i>AnnealingMethod</i>)
/// @param [in] chains The number of chain steps
/// @param [in] tolerancen Tolerance value to quit annealing
altar::bayesian::Annealer::
Annealer(
         problem_t & problem,
         scheduler_t & scheduler,
         sampler_t & sampler,
         method_t & worker,
         size_t chains, TYPE tolerance
         ) :
    _problem(problem), _scheduler(scheduler), _sampler(sampler), _worker(worker),
    _chains(chains),
    _tolerance(tolerance)
{
}


// end of file
