// -*- C++ -*-
// 
// michael a.g. aïvázis
// california institute of technology
// (c) 2010-2013 all rights reserved
// 


// for the build system
#include <portinfo>
#include <iomanip>
#include <iostream>

#include <altar/problem/Problem.h>

// get my declarations
#include "MPIAnnealing.h"
// externals
#include <altar/debug/print.h>


/// @par Main functionality
/// Update model data (eg gf, d) for Cp implementation
/// @param [in, out] annealer Host of the problem entity where priors, models and data reside
/// @return True if Cp implemented, false if Cp is not implemented
bool
altar::bayesian::MPIAnnealing::
updateModelData(annealer_t & annealer)
{
        /* problem_t & problem = static_cast<problem_t &>(annealer.problem()); */
        annealer_t::problem_t & problem = annealer.problem();
        return(problem.updateModelData(annealer, _worker.state().samples()));
}


/// @par Main functionality
/// Start the annealing process from scratch<br>
/// @param [in, out] annealer Host of the problem entity where priors, models and data reside
/// @param [in] fprior File name of the to-be-read prior
altar::bayesian::AnnealingMethod &
altar::bayesian::MPIAnnealing::
start(annealer_t & annealer, const std::string & fprior)
{
    // chain up
    AnnealingMethod::start(annealer, fprior);
    // notify my worker
    _worker.start(annealer, fprior);
    // collect the global state
    _collect();
    // all done
    return *this;
}

/// @par Main functionality
/// Start the annealing process from a previous run<br>
/// @param [in, out] annealer Host of the problem entity where priors, models and data reside
/// @param [in] iteration The iteration number to be restarted
altar::bayesian::AnnealingMethod &
altar::bayesian::MPIAnnealing::
restart(annealer_t & annealer, const size_t iteration)
{
    // chain up
    AnnealingMethod::restart(annealer, iteration);
    // notify my worker
    _worker.restart(annealer, iteration);
    // collect the global state
    _collect();
    // all done
    return *this;
}

/// @par Main functionality
/// Push my state forward along the cooling schedule
/// @param [in, out] annealer Host of the problem entity where priors, models and data reside
altar::bayesian::AnnealingMethod &
altar::bayesian::MPIAnnealing::
cool(annealer_t & annealer)
{
    // only the master task does anything useful
    if (!_master) return *this;
    // chain up
    return AnnealingMethod::cool(annealer);
}


/// @par Main functionality
/// Re-sample the posterior distribution
/// @param [in, out] annealer Host of the problem entity where priors, models and data reside
altar::bayesian::AnnealingMethod::stats_t
altar::bayesian::MPIAnnealing::
resample(annealer_t & annealer)
{
    // distribute my state
    _partition();
    // delegate the gathering of statistics to my worker
    stats_t stats = _worker.resample(annealer);
    // collect my state
    _collect();
    // and return the statistics
    return stats;
}

/// @par Main functionality
/// Analyze the acceptance statistics and take the problem state to the end of the current annealing step
/// @param [in, out] annealer Host of the problem entity where priors, models and data reside
/// @param [in, out] stats Alias for (acceptedSamples, rejectedSamples)
altar::bayesian::AnnealingMethod &
altar::bayesian::MPIAnnealing::
equilibrate(annealer_t & annealer, const annealer_t::stats_t & stats)
{
    // gather the statistics from every task
    stats_t globalStats = _gatherStatistics(stats);
    // braodcast from the root to all the nodes
    _bcastStatistics(globalStats);
    // if this is not the master task, we are done
//    if (!_master) return *this;
    // on the master task, chain up
    return AnnealingMethod::equilibrate(annealer, globalStats);
}

/// @par Main functionality
/// Notify me when annealing is finished
/// @param [in] annealer Host of the problem entity where priors, models and data reside
altar::bayesian::AnnealingMethod &
altar::bayesian::MPIAnnealing::
finish(annealer_t & annealer)
{
    // notify my worker
    _worker.finish(annealer);
    // chain up
    return AnnealingMethod::finish(annealer);
}

/// @par Main functionality
/// Status report
/// @param [in, out] annealer Host of the problem entity where priors, models and data reside
altar::bayesian::AnnealingMethod &
altar::bayesian::MPIAnnealing::
report(annealer_t & annealer)
{
    // only the master reports
    if (!_master) return *this;

    // chain up
    return AnnealingMethod::report(annealer);
}

/// @par Main functionality
/// Write the annealing progress to the "results" folder
/// @param [in] annealer Host of the problem entity where priors, models and data reside
/// @param [in, out] stats Alias for (acceptedSamples, rejectedSamples)
/// @param [in] restart The restart iteration number
altar::bayesian::AnnealingMethod &
altar::bayesian::MPIAnnealing::
write_out(annealer_t & annealer, const stats_t & stats, const bool flag_detail, const size_t restart)
{
    // only the master reports
    if (!_master) return *this;

    // chain up
    return AnnealingMethod::write_out(annealer, stats, flag_detail, restart);
}


/// @par Main functionality
/// destructor
altar::bayesian::MPIAnnealing::
~MPIAnnealing() 
{}


/// @par Main functionality
/// status report engine
/// @param [in, out] annealer Host of the problem entity where priors, models and data reside
altar::bayesian::AnnealingMethod &
altar::bayesian::MPIAnnealing::
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
    // show me contents
    altar::debug::print(_info, m, indent);

    // show me the prior
    state_t::vector_t * v = _state.prior();
    _info
        << "          prior: " << v->size << " samples" << pyre::journal::newline;
    // show me contents
    altar::debug::print(_info, v, indent);

    // show me the data
    v = _state.data();
    _info
        << "           data: " << v->size << " samples" << pyre::journal::newline;
    // show me contents
    altar::debug::print(_info, v, indent);

    // show me the posterior
    v = _state.posterior();
    _info
        << "      posterior: " << v->size << " samples" << pyre::journal::newline;
    // show me contents
    altar::debug::print(_info, v, indent);

    // show me {sigma}
    m = _state.sigma();
    _info
        << "          sigma: (" << m->size1 
        << " parameters) x (" << m->size2 << " parameters)" << pyre::journal::newline;
    // show me contents
    altar::debug::print(_info, m, indent);

    // all done
    return *this;
}

/// @par Main functionality
/// collect the states from MPI workers
void
altar::bayesian::MPIAnnealing::
_collect()
{
    // get the step information from my worker
    const state_t & localState = _worker.state();
    // collect the samples
    _gatherMatrix(_state.theta(), localState.theta());
    // collect the prior, data, and posterior
    _gatherVector(_state.prior(), localState.prior());
    _gatherVector(_state.data(), localState.data());
    _gatherVector(_state.posterior(), localState.posterior());

    // broadcast the global state
    MPI_Bcast(_state.theta()->data, _state.theta()->size1 * _state.theta()->size2, MPI_TYPE, 0, _communicator.handle());
    // copy the parameter covariance
//    gsl_matrix_TYPE_memcpy(_state.sigma(), localState.sigma());

    // all done
    return;
}

/// @par Main functionality
/// partition the states to MPI workers
void
altar::bayesian::MPIAnnealing::
_partition()
{
    // get the step information from my worker
    state_t & localState = _worker.state();
    
    // broadcast the temperature
    TYPE beta = _master ? _state.beta() : 0;
    MPI_Bcast(&beta, 1, MPI_TYPE, 0, _communicator.handle());
    localState.beta(beta);
    _state.beta(beta);

    // scatter the sample matrix
    _scatterMatrix(localState.theta(), _state.theta());
    // scatter prior, data, and posterior
    _scatterVector(localState.prior(), _state.prior());
    _scatterVector(localState.data(), _state.data());
    _scatterVector(localState.posterior(), _state.posterior());

    // copy the parameter covariance
    _bcastMatrix(localState.sigma(), _state.sigma());

    // all done
    return;
}

/// @par Main functionality
/// statistics gathering
void
altar::bayesian::MPIAnnealing::
_gatherMatrix(matrix_t * destination, matrix_t * source)
{
    // the place to deposit the data
    TYPE * data = _master? destination->data : 0;
    // the length of each contribution
    int size = source->size1 * source->size2;

    // gather
    MPI_Gather(
              source->data, size, MPI_TYPE, // send buffer
              data, size, MPI_TYPE, // receive buffer
              0, _communicator.handle() // address
              );

    // all done
    return;
}

/// @par Main functionality
/// statistics broadcasting
void
altar::bayesian::MPIAnnealing::
_bcastMatrix(matrix_t * destination, matrix_t * source)
{

    // at the master task
    if (_master) {
        // copy the source to the destination
        gsl_matrix_TYPE_memcpy(destination, source);
    }
    // compute the size
    size_t size = source->size1 * source->size2;
    // broadcast
    MPI_Bcast(destination->data, size, MPI_TYPE, 0, _communicator.handle());
    // all done
    return;
}

/// @par Main functionality
/// scatter gsl_matrix via MPI_Scatter
void
altar::bayesian::MPIAnnealing::
_scatterMatrix(matrix_t * destination, matrix_t * source)
{
    // the place to deposit the data
    TYPE * data = _master? source->data : 0;
    // the length of each contribution
    int size = destination->size1 * destination->size2;

    // gather
    MPI_Scatter(
              data, size, MPI_TYPE, // send buffer
              destination->data, size, MPI_TYPE, // receive buffer
              0, _communicator.handle() // address
              );

    // all done
    return;
}

/// @par Main functionality
/// gether gsl_vector via MPI_Gather
void
altar::bayesian::MPIAnnealing::
_gatherVector(vector_t * destination, vector_t * source)
{
    // the place to deposit the data
    TYPE * data = _master? destination->data : 0;
    // the length of each contribution
    int size = source->size;

    // gather
    MPI_Gather(
              source->data, size, MPI_TYPE, // send buffer
              data, size, MPI_TYPE, // receive buffer
              0, _communicator.handle() // address
              );

    // all done
    return;
}

/// @par Main functionality
/// scatter gsl_vector via MPI_Scatter
void
altar::bayesian::MPIAnnealing::
_scatterVector(vector_t * destination, vector_t * source)
{
    // the place to deposit the data
    TYPE * data = _master? source->data : 0;
    // the length of each contribution
    int size = destination->size;

    // gather
    MPI_Scatter(
              data, size, MPI_TYPE, // send buffer
              destination->data, size, MPI_TYPE, // receive buffer
              0, _communicator.handle() // address
              );

    // all done
    return;
}

/// @par Main functionality
/// statistics gathering
altar::bayesian::MPIAnnealing::stats_t
altar::bayesian::MPIAnnealing::
_gatherStatistics(const stats_t & stats)
{
    // unwrap
    TYPE accepted = stats.first;
    TYPE rejected = stats.second;

    // the resting place for the global stats
    TYPE totalAccepted, totalRejected;

    // collect
    MPI_Reduce(&accepted, &totalAccepted, 1, MPI_TYPE, MPI_SUM, 0, _communicator.handle());
    MPI_Reduce(&rejected, &totalRejected, 1, MPI_TYPE, MPI_SUM, 0, _communicator.handle());
    
    // all done
    return stats_t(totalAccepted,totalRejected);
}

/// @par Main functionality
/// statistics broadcasting
void
altar::bayesian::MPIAnnealing::
_bcastStatistics(stats_t & stats)
{
    // broadcast 
    MPI_Bcast(&stats.first, 1, MPI_TYPE, 0, _communicator.handle());
    MPI_Bcast(&stats.second, 1, MPI_TYPE, 0, _communicator.handle());
}


// end of file
