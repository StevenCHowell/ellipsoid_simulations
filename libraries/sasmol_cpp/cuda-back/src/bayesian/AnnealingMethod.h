// -*- C++ -*-
// -*- coding: utf-8 -*-
//
// michael a.g. aïvázis
// california institute of technology
// (c) 2010-2013  all rights reserved
//

// code guard
#if !defined(altar_bayesian_AnnealingMethod_h)
#define altar_bayesian_AnnealingMethod_h

// macros
#include <altar/utils/common.h>

// place everything in the local namespace
namespace altar {
    namespace problem {
        class CoolingStep;
    }
    namespace bayesian {

        // forward declarations
        class AnnealingMethod;

    } // of namespace bayesian
} // of namespace altar

// externals
#include <string>
#include <iostream>
#include <fstream>
#include <dirent.h>
#include <sys/stat.h>
#include <iomanip>
#include <pyre/journal.h>
#include <vector>
#ifdef USE_DOUBLE 
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#else
#include <gsl/gsl_vector_float.h>
#include <gsl/gsl_matrix_float.h>
#endif

// my dependencies
#include "Annealer.h"
#include <altar/problem/CoolingStep.h>


/// @brief Base class for the various annealing methods
///
/// @par Primary contents
/// CoolingStep: the place where @b M, <i>LLK</i>'s and @b Cs resides<br>
///
/// @par Main functionalities
/// start the annealing process from scratch<br>
/// restart the annealing process from a previous run<br>
/// re-sample the posterior distribution<br>
/// update model data for Cp implementation<br>
///
/// @note
/// This is the home of CoolingStep: the place where @b M, <i>LLK</i>'s and @b Cs resides<br>
/// <i> Note to developers:</i><br>
/// none
/// 

class altar::bayesian::AnnealingMethod
{
    // types
public:
    typedef Annealer annealer_t;///< alias for Annealer
    typedef Annealer::stats_t stats_t;///< alias for (acceptedSamples, rejectedSamples)
    typedef Annealer::sampler_t sampler_t;///< alias for Metropolis
    typedef Annealer::scheduler_t scheduler_t;///< alias for COV
    typedef altar::problem::CoolingStep state_t;///< alias for CoolingStep
    typedef AnnealingMethod method_t;///< alias for AnnealingMethod

    // accessors
public:
    inline TYPE beta() const;///< get the annealing temperature
    inline size_t iteration() const;///< get the iteration count of the present Annealing step 

    inline state_t & state();///< get reference of CoolingStep
    inline const state_t & state() const;///< get reference of CoolingStep as a constant

    // interface
public:
    /// update model data for Cp implementation
    virtual bool updateModelData(annealer_t &);
    /// start the annealing process from scratch
    virtual AnnealingMethod & start(annealer_t &, const std::string &);
    /// restart the annealing process from a previous run
    virtual AnnealingMethod & restart(annealer_t &, const size_t);
    /// notify me when annealing is finished
    virtual AnnealingMethod & finish(annealer_t &);

    // the annealing steps
    /// push my state forward along the cooling schedule
    virtual AnnealingMethod & cool(annealer_t &);
    /// re-sample the posterior distribution
    virtual stats_t resample(annealer_t &);
    /// analyze the acceptance statistics and take the problem state to the end of the current annealing step
    virtual AnnealingMethod & equilibrate(annealer_t &, const stats_t &);

    /// status report
    virtual AnnealingMethod & report(annealer_t &);
    /// write the annealing progress to the "results" folder
    virtual AnnealingMethod & write_out(annealer_t &, const stats_t & stats, const bool flag_detail, const size_t restart=0);

    // meta-methods
public:
    /// destructor
    virtual ~AnnealingMethod();
    /// constructor
    inline AnnealingMethod(std::string name, size_t workers, size_t samples, size_t parameters);
    /// iteration counter increment
    void iterate();

    // implementation details
protected:
    /// status report engine
    virtual AnnealingMethod & _report(annealer_t &);
private:
    /// read beta and scaling factor from the statistics file
    virtual AnnealingMethod & _readBetastatistics(std::string file, const size_t iteration, TYPE & beta, TYPE & scaling);

    // data
protected:
    std::string _name;///< name of the annealing method (eg distributed annealing, cuda annealing)
    size_t _workers;///< number of slaves that the present AnnealingMethod owns

    size_t _iteration;///< the iteration count of the present Annealing step
    state_t _state;///< CoolingStep: the place where @b M, <i>LLK</i>'s and @b Cs resides

    pyre::journal::info_t _info;///< debug stuff
    pyre::journal::debug_t _debug;///< debug stuff

    // disallow
private:
    /// copy constructor disallowed
    inline AnnealingMethod(const AnnealingMethod &);
    /// assign constructor disallowed
    inline const AnnealingMethod & operator=(const AnnealingMethod &);
};

// get the inline definitions
#define altar_bayesian_AnnealingMethod_icc
#include "AnnealingMethod.icc"
#undef altar_bayesian_AnnealingMethod_icc

# endif
// end of file
