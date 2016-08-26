// -*- C++ -*-
// -*- coding: utf-8 -*-
//
// michael a.g. aïvázis
// california institute of technology
// (c) 2010-2013  all rights reserved
//

// code guard
#if !defined(altar_bayesian_SequentialAnnealing_h)
#define altar_bayesian_SequentialAnnealing_h

// macros
#include <altar/utils/common.h>

// place everything in the local namespace
namespace altar {
    namespace problem {
        class CoolingStep;
    }
    namespace bayesian {

        // forward declarations
        class SequentialAnnealing;

    } // of namespace bayesian
} // of namespace altar

// my superclass
#include "AnnealingMethod.h"
// my dependencies
#include "Annealer.h"

// external
#include <string>

/// @brief Sequential annealing method
///
/// @par Primary contents
/// none
///
/// @par Main functionalities
/// start the annealing process from scratch<br>
/// restart the annealing process from a previous run<br>
///
/// @note
/// none
/// 

class altar::bayesian::SequentialAnnealing :
    public altar::bayesian::AnnealingMethod
{
    // types
public:
    typedef Metropolis sampler_t;///< alias for Metropolis
    typedef altar::problem::CoolingStep state_t;///< alias for CoolingStep
    typedef altar::problem::Problem problem_t;///< alias for Problem

    // interface
public:
    /// start the annealing process from scratch
    virtual AnnealingMethod & start(annealer_t &, const std::string &);
    /// restart the annealing process from a previous run
    virtual AnnealingMethod & restart(annealer_t &, const size_t);

    // meta-methods
public:
	/// constructor
    virtual ~SequentialAnnealing();
	/// destructor
    inline SequentialAnnealing(size_t samples, size_t parameters, unsigned short rank=0);

	// private member
private:
	const unsigned short _rank;///< the MPI rank <br>@note why should Sequential annealing be aware of the MPI rank?

    // implementation details
protected:
    /// status report engine
    virtual AnnealingMethod & _report(annealer_t &);

    // disallow
private:
	/// copy constructor disallowed
    inline SequentialAnnealing(const SequentialAnnealing &);
	/// assign constructor disallowed
    inline const SequentialAnnealing & operator=(const SequentialAnnealing &);
};

// get the inline definitions
#define altar_bayesian_SequentialAnnealing_icc
#include "SequentialAnnealing.icc"
#undef altar_bayesian_SequentialAnnealing_icc

# endif
// end of file
