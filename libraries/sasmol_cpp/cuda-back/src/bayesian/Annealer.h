// -*- C++ -*-
// -*- coding: utf-8 -*-
//
// michael a.g. aïvázis
// Hailiang Zhang
// california institute of technology
// (c) 2010-2013  all rights reserved
//

// code guard
#if !defined(altar_bayesian_Annealer_h)
#define altar_bayesian_Annealer_h

// macros
#include <altar/utils/common.h>

// place everything in the local namespace
namespace altar {
    namespace bayesian {

        // forward declarations
        class Annealer;
        class AnnealingMethod;
        class COV;

    } // of namespace bayesian
    namespace problem {
        class Problem;
    } // of namespace problem
} // of namespace altar

// external
#include <string>

// my dependencies
//#include <altar/problem/Problem.h>
//#include <altar/problem/CoolingStep.h>
#include "Metropolis.h"

// utils
#include <altar/utils/util.h>

/*! \mainpage
  * @par Introduction
  * CATMIP stands fro "Cascading Adaptive Tempered Metropolis In Parallel".<br>
  * This technique, based on Transitional Markov chain Monte Carlo, makes it possible to sample distributions in many hundreds of dimensions, if the forward model is fast, or to sample computationally expensive forward models in smaller numbers of dimensions.<br>
  * The design of the algorithm is independent of the model being sampled, so CATMIP can be applied to many areas of research.<br>
  * @par Reference
  * S. E. Minson, M. Simons and J. L. Beck., 2013. Bayesian inversion for finite fault earthquake source models I—theory and algorithm. <i>Geophysical Journal International</i>, <b>194</b>(3), 1701-1726
  * @par Programming Model
  * @image html ProgrammingModel.png
  */
///
/// @brief The driver of altar
///
/// @par Primary contents
/// References of Problem<br>
/// References of COV<br>
/// References of Metropolis<br>
/// References of AnnealingMethod<br>
///
/// @par Main functionalities
/// sampling the posterior distribution of my problem<br>
/// status report
///
/// @note
/// <i> Note to developers:</i><br>
/// none
/// 

// declaration
class altar::bayesian::Annealer
{
    // types
public:
    typedef AnnealingMethod method_t;///< alias for AnnealingMethod
    typedef altar::problem::Problem problem_t;///< alias for Problem
    typedef COV scheduler_t;///< alias for COV
    typedef Metropolis sampler_t;///< alias for Metropolis
    typedef Metropolis::stats_t stats_t;///< alias for (acceptedSamples, rejectedSamples)

    // accessors
public:
    inline size_t samples() const;///< get the number of samples
    /// @todo change to an inline function
    size_t parameters() const; ///< get the number of parameters

    inline problem_t & problem() const;///< get Problem reference
    inline method_t & worker() const;///< get AnnealingMethod reference
    inline sampler_t & sampler() const;///< get Metropolis reference
    inline scheduler_t & scheduler() const;///< get COV reference

    // interface
public:
    /// sampling the posterior distribution of my problem
    virtual int posterior(const size_t restart=0, const size_t wfreq=0, const std::string & fprior=std::string());

    /// status report
    virtual Annealer & report();


    // meta-methods
public:
    /// destructor
    virtual ~Annealer();
    /// constructor
    Annealer(
             // parts
             problem_t &, scheduler_t &, sampler_t &, method_t &,
             // settings
             size_t chains=2, TYPE tolerance=0.005);

    // data
protected:
    // parts
    problem_t & _problem;///< Reference of Problem
    scheduler_t & _scheduler;///< Reference of COV
    sampler_t & _sampler;///< Reference of Metropolis
    method_t & _worker;///< Reference of Annealing Method

    // settings
    size_t _chains;///< number of chains
    TYPE _tolerance;///< tolerance value to quit annealing

    // disallow
private:
    /// copy constructor disallowed
    inline Annealer(const Annealer &);
    /// assign constructor disallowed
    inline const Annealer & operator=(const Annealer &);
};

// get the inline definitions
#define altar_bayesian_Annealer_icc
#include "Annealer.icc"
#undef altar_bayesian_Annealer_icc

# endif
// end of file
