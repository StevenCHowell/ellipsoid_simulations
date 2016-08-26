// -*- C++ -*-
// -*- coding: utf-8 -*-
//
// michael a.g. aïvázis
// Hailiang Zhang
// california institute of technology
// (c) 2010-2013  all rights reserved
//

// code guard
#if !defined(altar_bayesian_Metropolis_h)
#define altar_bayesian_Metropolis_h

// macros
#include <altar/utils/common.h>

// place everything in the local namespace
namespace altar {
	namespace problem {

		// forward declarations
        class CoolingStep;
        class Problem;

	} // of namespace problem
    namespace bayesian {

        // forward declarations
        class Metropolis;

    } // of namespace bayesian
} // of namespace altar

// externals
#include <utility>
#include <iostream>
#include <cmath>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_randist.h>
#ifdef USE_DOUBLE 
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#else
#include <gsl/gsl_vector_float.h>
#include <gsl/gsl_matrix_float.h>
#endif


/// @brief Metropolis random walk
///
/// @par Primary contents
/// the number of chain steps<br>
/// the scaling factor to be applied on @b Cs<br>
/// the empirical weights for acceptance/rejection<br>
/// the random number generator
///
/// @par Main functionalities
/// Random walk <br>
/// analyze the acceptance statistics and take the problem state to the end of the current annealing step
///
/// @note
/// <i> Note to developers:</i><br>
/// cudaMetropolis has its own curand generator and cublas handle, primarily used for random sample generation on GPU
///

// declaration
class altar::bayesian::Metropolis
{
    // types
public:
    typedef altar::problem::Problem problem_t;///< alias for Problem
    typedef altar::problem::CoolingStep state_t;///< alias for CoolingStep

    typedef gsl_rng rng_t;///< alias for gsl_rng
    typedef gsl_vector_TYPE vector_t;///< alias for gsl_vector
    typedef gsl_matrix_TYPE matrix_t;///< alias for gsl_matrix

    typedef std::pair<TYPE, TYPE> stats_t;///< alias for (acceptedSamples, rejectedSamples)

    // interface
public:
	/// Random walk
    virtual stats_t sample(state_t &, const problem_t &) const;
	/// analyze the acceptance statistics and take the problem state to the end of the current annealing step
    virtual void equilibrate(const stats_t & stats);

    // meta-methods
public:
	/// destructor
    virtual ~Metropolis();
	/// constructor
    inline Metropolis(
                      rng_t *, // the random number generator
                      size_t steps = 20, // the length of my Markov chains
                      TYPE scaling = .1, // factor by which the covariance matrix is scaled
                      TYPE acceptanceWeight = 8, // the relative weight of accepted samples
                      TYPE rejectionWeight = 1 // the relative weight of rejected samples
                      );
                      
    // implementation details
protected: // Hailiang changed from private
	/// set up the random walk by giving each parameter a normally distributed random offset
    void displace(const state_t & state, const matrix_t * chol) const;
	/// replace rejected samples in {candidate} with copies from {current}
    void filter(const state_t & current, const state_t & candidate, const vector_t * rejects) const;

    // data
 	/// @note Hailiang changed from private to protected for <c>_steps, _scaling, _acceptanceWeight, _rejectionWeight, _rng</c>
protected:
    size_t _steps;///< the number of chain steps
    TYPE _scaling;///< the scaling factor to be applied on @b Cs
    TYPE _acceptanceWeight, _rejectionWeight;///< (_acceptanceWeight, _rejectionWeight) the empirical weights for acceptance/rejection
    rng_t * const _rng; ///< the random number generator; once set, it's fixed

    // disallow
private:
	/// copy constructor disallowed
    inline Metropolis(const Metropolis &);
	/// assign constructor disallowed
    inline const Metropolis & operator=(const Metropolis &);

friend class AnnealingMethod;///< @note <i>Hailiang added</i><br>allow write_out method to access its members 
};

// get the inline definitions
#define altar_bayesian_Metropolis_icc
#include "Metropolis.icc"
#undef altar_bayesian_Metropolis_icc

# endif
// end of file
