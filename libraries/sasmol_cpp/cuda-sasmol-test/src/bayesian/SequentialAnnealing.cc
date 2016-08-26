// -*- C++ -*-
// 
// michael a.g. aïvázis
// california institute of technology
// (c) 2010-2013 all rights reserved
// 


// for the build system
#include <portinfo>

#include <altar/problem/Problem.h>

// get my declarations
#include "SequentialAnnealing.h"

#include "rank.h"


/// @par Main functionality
/// Start the annealing process from scratch<br>
///- sample is initialized<br>
///- likelihoods are calculated
/// @param [in, out] annealer Host of the problem entity where priors, models and data reside
/// @param [in] fprior File name of the to-be-read prior
altar::bayesian::AnnealingMethod &
altar::bayesian::SequentialAnnealing::
start(annealer_t & annealer, const std::string & fprior)
{
    // get the problem
    annealer_t::problem_t & problem = annealer.problem();

    // initialize the prior/model-ensemble (prior/data-LLK and rand-seeds
    problem.initialEnsemble(_state.samples(), altar::utils::random_seed());

    // ask it to initialize/read my sample set
    if (!fprior.length()) problem.initialSample(_state, _rank);
    else problem.readSample(_state, fprior, _rank);

    // and compute the log likelihoods
    problem.likelihoods(_state);

    // print-outs
    int rank = my_rank();
    if (rank==0)
    {
        problem.show(_state.samples());
        _state.show();
    }

    // all done
    return AnnealingMethod::start(annealer, fprior);
}


/// @par Main functionality
/// Start the annealing process from a previous run<br>
/// @param [in, out] annealer Host of the problem entity where priors, models and data reside
/// @param [in] iteration The iteration number to be restarted
/// @todo
/// this not not done yet!
altar::bayesian::AnnealingMethod &
altar::bayesian::SequentialAnnealing::
restart(annealer_t & annealer, const size_t iteration)
{
    // get the rank
    int rank = my_rank();

    // get the problem
    problem_t & problem = static_cast<problem_t &>(annealer.problem());

    // initialize the prior/model-ensemble (prior/data-LLK and rand-seeds
    problem.initialEnsemble(_state.samples(), altar::utils::random_seed());

    // ask it to read the restart samples, llks and others
    _state.initialRestart(iteration, rank);

    // print-outs
    if (rank==0)
    {
        problem.show(_state.samples());
        _state.show();
    }

    // all done
    return AnnealingMethod::restart(annealer, iteration);
}

/// @par Main functionality
/// Destructor
altar::bayesian::SequentialAnnealing::
~SequentialAnnealing() 
{}


/// @par Main functionality
/// Report progress
/// @param [in] annealer Host of the problem entity where priors, models and data reside
altar::bayesian::AnnealingMethod &
altar::bayesian::SequentialAnnealing::
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
