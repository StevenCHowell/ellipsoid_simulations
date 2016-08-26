// -*- C++ -*-
// 
// michael a.g. aïvázis
// Hailiang Zhang
// california institute of technology
// (c) 2010-2013 all rights reserved
// 


// for the build system
#include <portinfo>

// get my declaration
#include "Metropolis.h"

// my dependencies
#include <altar/problem/Problem.h>
#include <altar/problem/CoolingStep.h>

#include "rank.h"

/// @par Main functionality
/// Random walk
/// @param [in, out] current The CPU state where @b M, <i>LLK</i>'s and @b Cs resides
/// @param [in] problem The problem entity (including priors, models, data)
altar::bayesian::Metropolis::stats_t
altar::bayesian::Metropolis::
sample(state_t & current, const problem_t & problem) const
{
    // some profiling
    struct timeval t1, t2, t1all, t2all;
    bool P_altar_profiling = false, P_altar_metropolis_detail_profiling = false;
    if (getenv("ALTAR_PROFILING") && strcmp(getenv("ALTAR_PROFILING"),"true")) P_altar_profiling = true;
    if (getenv("ALTAR_METROPOLIS_DETAIL_PROFILING") && strcmp(getenv("ALTAR_METROPOLIS_DETAIL_PROFILING"),"true")) P_altar_metropolis_detail_profiling = true;
    int rank = my_rank();

    // get the problem sizes
    const size_t samples = current.samples();
    const size_t parameters = current.parameters();
    // grab the covariance matrix
    matrix_t * sigma = current.sigma();
    // allocate room for my Cholesky decomposed covariance matrix
    matrix_t * sigma_chol = gsl_matrix_TYPE_alloc(parameters, parameters);
    // make a copy
    gsl_matrix_TYPE_memcpy(sigma_chol, sigma);
    // scale it
    gsl_matrix_TYPE_scale(sigma_chol, _scaling);
    // Cholesky decompose
#ifdef USE_DOUBLE
    gsl_linalg_cholesky_decomp(sigma_chol);
#else
    gsl_matrix * sigma_chol_double = altar::utils::gsl_matrix_double_from_float(sigma_chol);
    gsl_linalg_cholesky_decomp(sigma_chol_double);
    altar::utils::gsl_matrix_set_float_from_double(sigma_chol, sigma_chol_double);
    gsl_matrix_free(sigma_chol_double);
#endif

    // initialize the counts
    size_t acceptedSamples = 0;
    size_t rejectedSamples = 0;
    // allocate the rejection mask vector
    vector_t * rejects = gsl_vector_TYPE_alloc(samples);
    // allocate a vector to hold the difference between the two posterior likelihoods
    vector_t * diff = gsl_vector_TYPE_alloc(samples);
    // allocate a vector with random numbers for the Metropolis acceptance
    vector_t * dice = gsl_vector_TYPE_alloc(samples);
    
    if (P_altar_profiling) gettimeofday(&t1,NULL);
    // walk the chains
    for (size_t link=0; link<_steps; ++link) {
        // make a candidate step
        state_t candidate(samples, parameters, current.beta());

        // copy the sample set over to the candidate
        gsl_matrix_TYPE_memcpy(candidate.theta(), current.theta());
        // perform a random jump
        displace(candidate, sigma_chol);
        if (P_altar_metropolis_detail_profiling) altar::utils::print_profiling(rank,t1, t2, "Metropolis-displace");

        // the random displacement may have generated candidates that are outside the
        // support of the problem, so we must give the problem an opportunity to reject them
        // clear out the mask
        gsl_vector_TYPE_set_zero(rejects);
        // hand our candidate and the rejection mask to the problem
        problem.verify(candidate, rejects);
        if (P_altar_metropolis_detail_profiling) altar::utils::print_profiling(rank,t1, t2, "Metropolis-verify");
 
        // restore the consistency of the candidate sample set by replacing the rejected
        // samples with copies from the original
        filter(current, candidate, rejects);

        // ask the problem to compute the log likelihoods of the candidate
        problem.likelihoods(candidate);
        if (P_altar_metropolis_detail_profiling) altar::utils::print_profiling(rank,t1, t2, "Metropolis-likelihoods");

        // fill the vector with the posterior differences with the candidate posterior
        gsl_vector_TYPE_memcpy(diff, candidate.posterior());
        // subtract the current one
        gsl_vector_TYPE_sub(diff, current.posterior());

        // fill {dice} with uniformly distribute probabilities
        for (size_t sample=0; sample<samples; ++sample) {
            // get a random number
            TYPE face = gsl_ran_flat(_rng, 0, 1);
            // place it in {dice}
            gsl_vector_TYPE_set(dice, sample, face);
        }

        // accept/reject: go through all the samples
        for (size_t sample=0; sample < samples; ++sample) {
            // a candidate is rejected if the problem considered it invalid
            bool rejected = gsl_vector_TYPE_get(rejects, sample) != 0;
            // or it is less likely that the original and has not been saved by the dice
            bool unlikely = std::log(gsl_vector_TYPE_get(dice, sample)) > gsl_vector_TYPE_get(diff, sample);
            // so
            if (rejected || unlikely) {
                // increment the relevant counter
                rejectedSamples += 1;
                // and move on
                continue;
            }

            // otherwise, update the acceptance count
            acceptedSamples += 1;

            // access the new sample
            gsl_vector_TYPE_view newSample = gsl_matrix_TYPE_row(candidate.theta(), sample);
            // and its likelihoods
            TYPE newPrior = gsl_vector_TYPE_get(candidate.prior(), sample);
            TYPE newData = gsl_vector_TYPE_get(candidate.data(), sample);
            TYPE newPosterior = gsl_vector_TYPE_get(candidate.posterior(), sample);

            // copy the sample over to the current state
            gsl_matrix_TYPE_set_row(current.theta(), sample, &newSample.vector);
            // and the likelihoods
            gsl_vector_TYPE_set(current.prior(), sample, newPrior);
            gsl_vector_TYPE_set(current.data(), sample, newData);
            gsl_vector_TYPE_set(current.posterior(), sample, newPosterior);
        }
        if (P_altar_metropolis_detail_profiling) altar::utils::print_profiling(rank,t1, t2, "Metropolis-update");
    }

    // free the vector temporaries
    gsl_vector_TYPE_free(diff);
    gsl_vector_TYPE_free(dice);
    gsl_vector_TYPE_free(rejects);
    // and the Cholesky decomposed {sigma}
    gsl_matrix_TYPE_free(sigma_chol);

    if (P_altar_profiling) altar::utils::print_profiling(rank,t1all, t2all, "Metropolis");

    // and return
    return stats_t(acceptedSamples, rejectedSamples);
}

/// @par Main functionality
/// Update my problem statistics based on the results of waling my Markov chains
/// @note
/// IMPLEMENTATION NOTE: this is factored out to accommodate the parallel versions that need to
/// collect the global state from each Markov chain walker before adjusting the covariance
/// matrix scaling
void
altar::bayesian::Metropolis::
equilibrate(const stats_t & stats)
{
    // alias my data
    TYPE aw = _acceptanceWeight;
    TYPE rw = _rejectionWeight;

    // unpack the statistics
    TYPE accepted = std::get<0>(stats);
    TYPE rejected = std::get<1>(stats);


    // compute the acceptance ratio
    TYPE acceptance =  accepted / (accepted+rejected);

    // compute the fudge factor
    TYPE kc = (aw*acceptance + rw)/(aw+rw);

    // don't let it get too small
    if (kc < .1) {
        kc = .1;
    // or too big
    } else if (kc > 1) {
        kc = 1;
    }

    // store it
    //_scaling = kc;
    _scaling = kc*kc;//Hailiang changed it


    // all done
    return;
}

/// @par Main functionality
/// destructor
altar::bayesian::Metropolis::
~Metropolis() 
{}


/// @par Main functionality
/// Set up the random walk by giving each parameter a normally distributed random offset
/// @param [in, out] candidate The CPU candidate state where @b M, <i>LLK</i>'s and @b Cs resides
/// @param [in] chol The sample covariance matrix after cholosky decomposition
void
altar::bayesian::Metropolis::
displace(const state_t & candidate, const matrix_t * chol) const
{
    // extract the problem sizes
    const size_t samples = candidate.samples();
    const size_t parameters = candidate.parameters();

    // build a set of normally distributed random vectors
    matrix_t * delta = gsl_matrix_TYPE_alloc(samples, parameters);
    // fill it
    for (size_t i=0; i<samples; ++i) {
        for (size_t j=0; j<parameters; ++j) {
            gsl_matrix_TYPE_set(delta, i, j, gsl_ran_ugaussian(_rng));
        }
    }

    // post-multiply by the transpose of the Cholesky decomposed sigma, which is stored in the
    // upper triangle of {chol}
    gsl_blas_TYPEtrmm(CblasRight, CblasUpper, CblasNoTrans, CblasNonUnit, 1, chol, delta);

    // increment the sample by delta
    gsl_matrix_TYPE_add(candidate.theta(), delta);

    // free the temporaries
    gsl_matrix_TYPE_free(delta);

    // all done
    return;
}

/// @par Main functionality
/// Replace rejected samples in <i>candidate</i> with copies from <i>current</i>
/// @param [in, out] current The CPU current state where @b M, <i>LLK</i>'s and @b Cs resides
/// @param [in] candidate The CPU candidate state where @b M, <i>LLK</i>'s and @b Cs resides
/// @param [in] rejects Pointer to the gsl_vector of the sample rejection flags
void
altar::bayesian::Metropolis::
filter(const state_t & current, const state_t & candidate, const vector_t * rejects) const
{
    // get the problem size
    const size_t samples = current.samples();

    // access the sample sets
    matrix_t * theta = current.theta(); // the original sample set
    matrix_t * prime = candidate.theta(); // the displace sample set

    // go through the mask
    for (size_t sample=0; sample<samples; ++sample) {
        // if this sample is rejected
        if (gsl_vector_TYPE_get(rejects, sample) != 0) {
            // access the corresponding row in {current}
            gsl_vector_TYPE_view original = gsl_matrix_TYPE_row(theta, sample);
            // copy it over to the {candidate}
            gsl_matrix_TYPE_set_row(prime, sample, &original.vector);
        }
    }
    // all done
    return;
}


// end of file
