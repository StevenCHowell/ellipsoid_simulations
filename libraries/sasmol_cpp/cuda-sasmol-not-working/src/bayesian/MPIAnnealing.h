// -*- C++ -*-
// -*- coding: utf-8 -*-
//
// michael a.g. aïvázis
// california institute of technology
// (c) 2010-2013  all rights reserved
//

// code guard
#if !defined(altar_bayesian_MPIAnnealing_h)
#define altar_bayesian_MPIAnnealing_h

// macros
#include <altar/utils/common.h>

// externals
#include <pyre/mpi.h>

// place everything in the local namespace
namespace altar {
    namespace bayesian {

        // forward declarations
        class MPIAnnealing;

    } // of namespace bayesian
} // of namespace altar

// my superclass
#include "AnnealingMethod.h"
// my dependencies
#include "Annealer.h"
#include <altar/problem/CoolingStep.h>

// external
#include <string>

/// @brief MPI annealing method
///
/// @par Primary contents
/// CoolingStep: the place where @b M, <i>LLK</i>'s and @b Cs resides<br>
/// reference of MPI worker
///
/// @par Main functionalities
/// start the annealing process from scratch<br>
/// restart the annealing process from a previous run<br>
/// re-sample the posterior distribution<br>
/// update model data for Cp implementation<br>
/// collect/partition the states from/to MPI workers
///
/// @note
/// This is the home of CoolingStep: the place where @b M, <i>LLK</i>'s and @b Cs resides<br>
/// <i> Note to developers:</i><br>
/// none
/// 

class altar::bayesian::MPIAnnealing :
    public altar::bayesian::AnnealingMethod
{
    // types
public:
    typedef altar::problem::CoolingStep::matrix_t matrix_t;///< alias for gsl_matrix
    typedef altar::problem::CoolingStep::vector_t vector_t;///< alias for gsl_vector
    typedef pyre::mpi::communicator_t communicator_t;///< alias for MPI communicator

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
    virtual AnnealingMethod & write_out(annealer_t &, const stats_t &, const bool flag_detail, const size_t restart=0);

    // getter
public:
    inline virtual communicator_t & communicator_const() const;///< get the communicator

    // meta-methods
public:
    /// destructor
    virtual ~MPIAnnealing();
    /// constructor
    inline MPIAnnealing(communicator_t & communicator,
                        size_t samples, size_t parameters, method_t & worker);

    // implementation details
protected:
    /// status report engine
    virtual AnnealingMethod & _report(annealer_t &);

    /// collect the states from MPI workers
    void _collect();
    /// partition the states to MPI workers
    void _partition();

    /// statistics gathering
    stats_t _gatherStatistics(const stats_t &);
    /// statistics broadcasting
    void _bcastStatistics(stats_t &);

    // mpi+gsl workhorses
    /// gether gsl_vector via MPI_Gather
    void _gatherVector(vector_t * destination, vector_t * source);
    /// scatter gsl_vector via MPI_Scatter
    void _scatterVector(vector_t * destination, vector_t * source);
    /// broadcast gsl_matrix via MPI_Bcast
    void _bcastMatrix(matrix_t * destination, matrix_t * source);
    /// gather gsl_matrix via MPI_Gather
    void _gatherMatrix(matrix_t * destination, matrix_t * source);
    /// scatter gsl_matrix via MPI_Scatter
    void _scatterMatrix(matrix_t * destination, matrix_t * source);

    // private data
private:
    bool _master; ///< a flag to mark the master task
    communicator_t & _communicator; ///< reference of my communicator
    method_t & _worker; ///< reference of my worker

    // disallow
private:
    /// copy constructor disallowed
    inline MPIAnnealing(const MPIAnnealing &);
    /// assign constructor disallowed
    inline const MPIAnnealing & operator=(const MPIAnnealing &);
};

// get the inline definitions
#define altar_bayesian_MPIAnnealing_icc
#include "MPIAnnealing.icc"
#undef altar_bayesian_MPIAnnealing_icc

# endif
// end of file
