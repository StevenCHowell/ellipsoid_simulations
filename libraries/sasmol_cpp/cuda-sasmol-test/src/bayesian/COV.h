// -*- C++ -*-
// -*- coding: utf-8 -*-
//
// michael a.g. aïvázis
// Hailiang Zhang
// california institute of technology
// (c) 2010-2013  all rights reserved
//

// code guard
#if !defined(altar_bayesian_COV_h)
#define altar_bayesian_COV_h

// macros
#include <altar/utils/common.h>

// externals
#include <cmath>
#include <iostream>
#include <iomanip>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_sys.h>
#include <gsl/gsl_min.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_histogram.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_sort_vector.h>
#include <gsl/gsl_permute_vector.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector_float.h>
#include <gsl/gsl_matrix_float.h>
#include <gsl/gsl_sort_vector_float.h>
#include <gsl/gsl_permute_vector_float.h>

// utils
#include <altar/bayesian/util.h>

// place everything in the local namespace
namespace altar {
    namespace problem {

        // forward declarations
        class CoolingStep;

    } // of namespace problem
    namespace bayesian {

        // forward declarations
        class COV;

    } // of namespace bayesian
} // of namespace altar

// workhorses
namespace cov {
    /// the COV calculator
    static double cov(double, void *);

    /// the auxiliary data structure for weight underflow control
    struct underflow_control_t
    {
        size_t index_begin; ///< the starting index of the scaled region
        double scale; ///< the scaling factor (>=1)
    };

    /// the parameter structure of COV calculator
    struct args {
        double dbeta; ///< the current proposal for dbeta
        double cov; ///< the current value of COV
        double metric; ///< the latest value of our objective function

        // inputs from the higher levels
        gsl_vector * w; ///< the vector of weights
        gsl_vector * llk; ///< the vector of data log-likelihoods
        double target; ///< the COV value we are aiming for; should be 1

        // auxiliary data for weight underflow control
        std::vector<underflow_control_t> wc;
    };

    /// scale the weights and save the scaling information
    static void scaleWeights(args &);
    /// return a gsl_vector that contains the fine weight information for each sample
    gsl_vector * fineWeights(const args &);
    /// return a gsl_vector that contains the unscaled coarse (sum) weight information for each sub_region
    gsl_vector * coarseWeights(const args &);
    /// get the sum of the scaled-weights beyond a given region (including this region)
    double sumWeights(const args & p, const int region_start=0);
    /// get the sum of the squared scaled-weights
    double sumWeightsSquared(const args & p);
}

/// @brief Cooling between beta steps
///
/// @par Primary contents
/// the annealing temperature<br>
/// the sample weight covariance<br>
/// the random number generator
///
/// @par Main functionalities
/// calculate the new annealing temperature<br>
/// compute the new annealing temperature<br>
/// update @b C<sub>s</sub><br>
/// reshuffle samples
///
/// @note
/// <i> Note to developers:</i><br>
/// @c rankAndShuffle function is not sorting the samples any more (due to MPI partition issue)
///
/// @bug
/// <i>dbeta_gsl</i> not working yet because there is no gsl minimizer for single precision
///

// declaration
class altar::bayesian::COV
{
    // types
public:
    typedef altar::problem::CoolingStep state_t;///< alias for CoolingStep

    typedef gsl_rng rng_t;///< alias for gsl_rng
    typedef gsl_vector_TYPE vector_t;///< alias for gsl_vector
    typedef gsl_matrix_TYPE matrix_t;///< alias for gsl_matrix

    // accessors
public:
    inline TYPE cov() const;///< get the covariance of the sample weights
    inline TYPE beta() const;///< get the annealing temperature

    // interface
public:
    virtual void update(state_t & state);///< <i>the primary method in COV:</i><br> compute the new annealing temperature, update @b C<sub>s</sub>, and reshuffle samples
    inline void beta(const TYPE);///< set the annealing temperature

    // lower level interface
public:
    virtual double dbeta_gsl(cov::args &);///< calculate the new annealing temperature by gsl minimizer
    virtual double dbeta_grid(cov::args &);///< calculate the new annealing temperature by iterative grid-based searching

    // meta-methods
public:
    /// destructor
    virtual ~COV();
    /// constructor
    inline COV(
               rng_t * rng,
               TYPE tolerance=.001, size_t maxIterations=1000, TYPE target=1.0
               );

    // implementation details
protected:
    void computeCovariance(state_t & state, const cov::args & covargs, const gsl_permutation * perm) const; ///< compute the new sample covariance matrix
    void computePmi(gsl_vector * Pmi, const cov::args & covargs, const gsl_permutation * perm) const; ///< compute the new sample covariance matrix
    void conditionCovariance(matrix_t * sigma) const;///< condition the covariance matrix
    void rankAndShuffle(state_t & state, const cov::args & covargs, const gsl_permutation * perm) const; ///< reshuffle samples

    // data
private:
    const TYPE _betaMin;///< the lower bound for dbeta optimization
    const TYPE _betaMax;///< the upper bound for dbeta optimization

    rng_t * _rng;///< random generator
    TYPE _tolerance;///< the tolerance value used in dbeta optimization that we can skip straight to beta = 1
    size_t _maxIterations;///< the maximum number of iterations for dbeta value optimization
    TYPE  _target;///< the COV value we are aiming for; should be 1
    TYPE _beta, _cov;///< (_beta, _cov) the annealing temperature and the sample weight covariance

    // disallow
private:
    /// copy constructor disallowed
    inline COV(const COV &);
    /// assign constructor disallowed
    inline const COV & operator=(const COV &);
};

// get the inline definitions
#define altar_bayesian_COV_icc
#include "COV.icc"
#undef altar_bayesian_COV_icc

# endif
// end of file
