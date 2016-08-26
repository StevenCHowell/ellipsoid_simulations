// -*- C++ -*-
// 
// michael a.g. aïvázis
// Hailiang Zhang
// california institute of technology
// (c) 2010-2013 all rights reserved
// 


// for the build system
#include <portinfo>

// for debugging
#include <cassert>


#include <pyre/journal.h>

// get my declarations
#include "COV.h"

// and my dependencies
#include <altar/problem/CoolingStep.h>

// my externels
#include <exception>
#include <limits>
#include <sys/stat.h>


/// @par Main functionality
/// Compute the new annealing temperature, update @b C<sub>s</sub>, and re-shuffle samples
/// @param [in, out] state The CPU state where @b M, <i>LLK</i>'s and @b Cs resides
/// @note
/// This is the primary method in COV
void
altar::bayesian::COV::
update(state_t & state)
{
    // get the problem size
    const size_t samples = state.samples();
    // copy the data likelihood vector
    gsl_vector * dataLLK = gsl_vector_alloc(samples);
#ifdef USE_DOUBLE
    gsl_vector_memcpy(dataLLK, state.data());
#else
    altar::utils::gsl_vector_set_double_from_float(dataLLK, state.data());
#endif
    // make a vector for the weights
    gsl_vector * w = gsl_vector_alloc(samples);

    // sort dataLLK in descending order and save the permutation
    gsl_permutation * p = gsl_permutation_alloc(samples);
    gsl_sort_vector_index(p, dataLLK);
    gsl_permute_vector(p, dataLLK);
    gsl_vector_reverse(dataLLK);
    gsl_permutation_reverse(p);

    // allocate space for the parameters
    cov::args covargs;
    // attach the two vectors
    covargs.w = w;
    covargs.llk = dataLLK;
    // initialize the COV target
    covargs.target = _target;

    // compute the temperature increment
    //dbeta_gsl(covargs);
    dbeta_grid(covargs);
    // save the new temperature
    state.beta(_beta);

    // compute the covariance
    computeCovariance(state, covargs, p);

    // rank and reorder the samples according to their likelihoods
    rankAndShuffle(state, covargs, p);
        
    // re-calcuate the posterior based on the new beta // Hailiang added
    vector_t * priorLLK = state.prior();
    vector_t * postLLK = state.posterior();
    gsl_blas_TYPEcopy(priorLLK, postLLK);
    gsl_blas_TYPEaxpy(state.beta(), state.data(), postLLK);

    // free the temporaries
    gsl_vector_free(dataLLK);
    gsl_vector_free(w);
    gsl_permutation_free(p);

    // all done
    return;
}

/// @par Main functionality
/// Calculate the new annealing temperature by gsl minimizer
/// @param [in] covargs Arguments for function minimization
/// @return The calculated beta value
double
altar::bayesian::COV::
dbeta_gsl(cov::args & covargs)
{
    // build my debugging channel
    pyre::journal::debug_t debug("altar.beta");

    // turn off the err_handler (Hailiang)
    gsl_error_handler_t * gsl_hdl = gsl_set_error_handler_off ();

    // consistency checks
    assert(_betaMin == 0);
    assert(_betaMax == 1);

    // the beta search region
    double beta_low = 0;
    double beta_high = _betaMax - _beta;
    double beta_guess = _betaMin + 5.0e-5;

    assert(beta_high >= beta_low);
    assert(beta_high >= beta_guess);
    assert(beta_guess >= beta_low);

    // search bounds and initial guess
    double f_beta_high = cov::cov(beta_high, &covargs);

    // check whether we can skip straight to beta = 1
    if (covargs.cov < _target || std::abs(covargs.cov-_target) < _tolerance) {
        debug 
            << pyre::journal::at(__HERE__)
            << " ** skipping to beta = " << _betaMax << " **" 
            << pyre::journal::endl;
        // save my state
        _beta = _betaMax;
        _cov = covargs.cov;
        // all done
        return beta_high;
    }

    double f_beta_low = cov::cov(beta_low, &covargs);
    // do this last so our first printout reflects the values at out guess
    double f_beta_guess = cov::cov(beta_guess, &covargs);

    // lie to the minimizer, if necessary
    // if (f_beta_low < f_beta_guess) f_beta_low = 1.01 * f_beta_guess;
    // if (f_beta_high < f_beta_guess) f_beta_high = 1.01 * f_beta_guess;

    // instantiate the minimizer
    gsl_min_fminimizer * minimizer = gsl_min_fminimizer_alloc(gsl_min_fminimizer_brent);
    // set up the minimizer call back
    gsl_function F;
    F.function = cov::cov;
    F.params = &covargs;
    // prime it
    gsl_min_fminimizer_set_with_values(
                                       minimizer, &F,
                                          beta_guess, f_beta_guess,
                                          beta_low, f_beta_low,
                                          beta_high, f_beta_high);
    // duplicate the catmip output, for now
    size_t iter = 0;
    if (debug) {
        std::cout
            << "    "
            << "Calculating dbeta using "
            << gsl_min_fminimizer_name(minimizer) << " method" << std::endl
            << "      target: " << _target << std::endl
            << "      tolerance: " << _tolerance << std::endl
            << "      max iterations: " << _maxIterations
            << std::endl;
        std::cout
            << "     "
            << std::setw(6) << " iter"
            " [" << std::setw(11) << "lower" 
            << ", " << std::setw(11) << "upper" << "]"
            << " " << std::setw(11) << "dbeta  " 
            << " " << std::setw(11) << "cov  " 
            << " " << std::setw(11) << "err  " 
            << " " << std::setw(11) << "f(dbeta)"
            << std::endl;
        std::cout.setf(std::ios::scientific, std::ios::floatfield);
        std::cout
            << "     "
            << std::setw(5) << iter
            << " [" << std::setw(11) << std::setprecision(4) << beta_low
            << ", " << std::setw(11) << beta_high << "] "
            << " " << std::setw(11) << covargs.dbeta
            << " " << std::setw(11) << covargs.cov
            << " " << std::setw(11) << covargs.cov - _target
            << " " << std::setw(11) << covargs.metric
            << std::endl;
    }

    // the GSL flag
    int status;
    // iterate, looking for the minimum
    do {
        iter++;
        status = gsl_min_fminimizer_iterate(minimizer);
        beta_low = gsl_min_fminimizer_x_lower(minimizer);
        beta_high = gsl_min_fminimizer_x_upper(minimizer);
        status = gsl_root_test_residual(covargs.cov-_target,  _tolerance);

        // print the result
        if (debug) {
            std::cout
                << "     "
                << std::setw(5) << iter
                << " [" << std::setw(11) << std::setprecision(4) << beta_low
                << ", " << std::setw(11) << beta_high << "] "
                << " " << std::setw(11) << covargs.dbeta
                << " " << std::setw(11) << covargs.cov
                << " " << std::setw(11) << covargs.cov - _target
                << " " << std::setw(11) << covargs.metric;
            
            if (status == GSL_SUCCESS) {
                std::cout << " (Converged)";
            }
            std::cout << std::endl;
        }

    } while (status == GSL_CONTINUE && iter < _maxIterations);

    // get the best guess at the minimum
    double dbeta = gsl_min_fminimizer_x_minimum(minimizer);
    // make sure that we are left with the COV and weights evaluated for this guess
    cov::cov(dbeta, &covargs);

    // free the minimizer
    gsl_min_fminimizer_free(minimizer);
    
    // adjust my state
    _cov = covargs.cov;
    _beta += dbeta;
    // return the beta update
    return dbeta;
}


/// @par Main functionality
/// Calculate the new annealing temperature by iterative grid-based searching
/// @param [in] covargs Arguments for function minimization
/// @return The calculated beta value
double
altar::bayesian::COV::
dbeta_grid(cov::args & covargs)
{
    // build my debugging channel
    pyre::journal::debug_t debug("altar.beta");

    // turn off the err_handler (Hailiang)
    gsl_error_handler_t * gsl_hdl = gsl_set_error_handler_off ();

    // consistency checks
    assert(_betaMin == 0);
    assert(_betaMax == 1);

    // the beta search region
    double beta_low = 0;
    double beta_high = _betaMax - _beta;
    double beta_guess = _betaMin + 5.0e-5;

    assert(beta_high >= beta_low);
    assert(beta_high >= beta_guess);
    assert(beta_guess >= beta_low);

    // search bounds and initial guess
    double f_beta_high = cov::cov(beta_high, &covargs);

    // check whether we can skip straight to beta = 1
    if (covargs.cov < _target || std::abs(covargs.cov-_target) < _tolerance) {
        debug 
            << pyre::journal::at(__HERE__)
            << " ** skipping to beta = " << _betaMax << " **" 
            << pyre::journal::endl;
        // save my state
        _beta = _betaMax;
        _cov = covargs.cov;
        // all done
        return beta_high;
    }

    double f_beta_low = cov::cov(beta_low, &covargs);
    // do this last so our first printout reflects the values at out guess
    double f_beta_guess = cov::cov(beta_guess, &covargs);

    // iterative grid searching
    //double beta_grid_tolerance=1.E-6;
    const int Nbeta=10;
    int nbeta = 0;
    double dbeta, beta_step;
    int Nloop = 0;
    if (debug) std::cout<<"dbeta minimization based on iterative grid searching:"<<std::endl;
    bool Qfind = false;
    do
    {
        ++Nloop;
        beta_step = (beta_high-beta_low)/Nbeta;
        int count=0;
        for (beta_guess=beta_low, nbeta=0; nbeta<=Nbeta; beta_guess+=beta_step, ++nbeta)
        {
            f_beta_guess = cov::cov(beta_guess, &covargs);
            if (std::abs(covargs.cov-_target) < _tolerance)
            {
                Qfind = true;
                break;
            }
            if (covargs.cov >= _target)
            {
               if (nbeta==0) Qfind = true;
               break;
            }
        }
        if (nbeta>Nbeta) beta_guess -= beta_step;
        if (debug) {
            std::cout 
                <<"      dbeta_low: "<<std::setw(11)<<std::setprecision(8)<<beta_low
                <<"      dbeta_high: "<<std::setw(11)<<beta_high
                <<"      dbeta_guess: "<<std::setw(11)<<beta_guess
                <<"      cov: "<<std::setw(15)<<covargs.cov
                << std::endl;
        }
        beta_high = beta_guess;
        beta_low = beta_guess-beta_step;
    }
    while (Nloop<=Nbeta && !Qfind);

    // assign dbeta
    dbeta = beta_guess;

    // make sure that we are left with the COV and weights evaluated for this guess
    cov::cov(dbeta, &covargs);
    
    // adjust my state
    _cov = covargs.cov;
    _beta += dbeta;
    // return the beta update
    return dbeta;
}


/// @par Main functionality
/// destructor
altar::bayesian::COV::
~COV() 
{}

/// @par Main functionality
/// Compute the new sample covariance matrix
/// @param [in, out] state The CPU state where @b M, <i>LLK</i>'s and @b Cs resides
/// @param [in] weights Pointer to the gsl_vector of the sample weights
void
altar::bayesian::COV::
computePmi(gsl_vector * Pmi, const cov::args & covargs, const gsl_permutation * perm) const
{
    // get the system size
    const size_t samples = Pmi->size;
    const size_t regions = covargs.wc.size();

    // get the weight and weight_control
    const gsl_vector * w = covargs.w;
    const std::vector<cov::underflow_control_t> & wc = covargs.wc;

    // local variables
    int region, sample, samples_local, index_begin;
    double pmi;
    double scale = 1.0;
    double sum_weightUnscaled_local, sum_weightScaled_beyond;

    // do the work
    // get the sum of scaled weights for all regions
    double sum_weightScaled_all = sumWeights(covargs);
    // sequentially loop over coarse_weight regions
    for (region=0; region<regions; ++region)
    {
        // get the starting index for this region
        index_begin = wc[region].index_begin;
        // get the scaling factor for this region w.r.t the 1st region
        scale *= wc[region].scale;
        // get the number of samples for this region
        if (region<regions-1) samples_local = wc[region+1].index_begin-index_begin;
        else samples_local = samples-index_begin;
        // get the sum of unscaled weights for this region only
        gsl_vector_const_view w_sub = gsl_vector_const_subvector(w, index_begin, samples_local);
        sum_weightUnscaled_local = gsl_blas_dasum(&w_sub.vector);
        // get the sum of scaled weights beyond this region (inclusive)
        sum_weightScaled_beyond = sumWeights(covargs, region);
        // set the scaled pmi for each sample in this region
        for (sample=0; sample<samples_local; ++sample)
        {
            pmi = gsl_vector_get(w, index_begin+sample)/sum_weightUnscaled_local;
            //pmi *= (sum_weightScaled_beyond/sum_weightScaled_all)*scale;
            pmi *= (sum_weightUnscaled_local*scale)/sum_weightScaled_all;
            gsl_vector_set(Pmi, index_begin+sample, pmi);
        }
    }
    
    // normalize Pmi
    double wsum = Pmi->size * gsl_stats_mean(Pmi->data, Pmi->stride, Pmi->size);
    gsl_vector_scale(Pmi, 1/wsum);

    // permutate Pmi
    gsl_permutation * perm_inv = gsl_permutation_alloc(samples);
    gsl_permutation_inverse(perm_inv, perm);
    gsl_permute_vector(perm_inv, Pmi);

    // free locals
    gsl_permutation_free(perm_inv);

    // all done
    return;
}


/// @par Main functionality
/// Compute the new sample covariance matrix
/// @param [in, out] state The CPU state where @b M, <i>LLK</i>'s and @b Cs resides
/// @param [in] weights Pointer to the gsl_vector of the sample weights
void
altar::bayesian::COV::
computeCovariance(state_t & state, const cov::args & covargs, const gsl_permutation * perm) const
{
    // unpack the problem sizes
    const size_t samples = state.samples();
    const size_t parameters = state.parameters();

    // get the scaled and re-odered pmi (weights)
    gsl_vector * weights_double = gsl_vector_alloc(samples);
    computePmi(weights_double, covargs, perm);    

    // make a vector for the weights
    vector_t * weights = gsl_vector_TYPE_alloc(samples);
#ifdef USE_DOUBLE
    gsl_vector_memcpy(weights, weights_double);
#else
    altar::utils::gsl_vector_set_float_from_double(weights, weights_double);
#endif

    // get the sample matrix
    matrix_t * theta = state.theta();
    // and the covariance
    matrix_t * sigma = state.sigma();
    
    // zero it out before we start accumulating values there
    gsl_matrix_TYPE_set_zero(sigma);

    // build a vector for the weighted mean of each parameter
    vector_t * thetaBar = gsl_vector_TYPE_alloc(parameters);
    // for each parameter
    for (size_t parameter=0; parameter<parameters; ++parameter) {
        // get column of theta that has the value of this parameter across all samples
        gsl_vector_TYPE_view column = gsl_matrix_TYPE_column(theta, parameter);
        // and treat it like a vector
        vector_t * values = &column.vector;
        // compute the mean
        TYPE mean = gsl_stats_TYPE_wmean(
                                      weights->data, weights->stride,
                                      values->data, values->stride, values->size
                                      );
        // set the corresponding value of theta bar
        gsl_vector_TYPE_set(thetaBar, parameter, mean);
    }

    // start filling out sigma
    for (size_t sample=0; sample<samples; ++sample) {
        // get the row with the sample
        gsl_vector_TYPE_view row = gsl_matrix_TYPE_row(theta, sample);
        // form {sigma += w[i] sample sample^T}
        gsl_blas_TYPEsyr(CblasLower, gsl_vector_TYPE_get(weights, sample), &row.vector, sigma);
    }
    // subtract {thetaBar . thetaBar^T}
    gsl_blas_TYPEsyr(CblasLower, -1, thetaBar, sigma);

    // fill the upper triangle
    for (size_t i=0; i<parameters; ++i) {
        for (size_t j=0; j<i; ++j) {
            gsl_matrix_TYPE_set(sigma, j,i, gsl_matrix_TYPE_get(sigma, i,j));
        }
    }

    // condition the covariance matrix
    //altar::utils::mprintf(sigma, "New Sigma before conditioning:", 15, 15);
    altar::bayesian::util::conditionMatrix(sigma);
    //altar::utils::mprintf(sigma, "Sigma after conditioning:", 15, 15);

    
    // free the temporaries
    gsl_vector_TYPE_free(thetaBar);
    gsl_vector_free(weights_double);
    gsl_vector_TYPE_free(weights);

    // all done
    return;
}

/// @par Main functionality
/// Compute the new sample covariance matrix
/// @param [in, out] sigma The sample covariance matrix
/// @note
/// This function is actually not implemented<br>
/// The actually matrix conditioning is implemented in "altar::bayesian::util::conditionMatrix"<br>
/// because this is a general matrix computation
void
altar::bayesian::COV::
conditionCovariance(matrix_t * sigma) const
{
    // all done
    return;
}

/// @par Main functionality
/// Reshuffle samples
/// @note
/// Ranking is not done despite the function name
/// This is because ranking (sorting) will produce bias for the subsequent MPI sample partitioning
/// @param [in, out] state The CPU state where @b M, <i>LLK</i>'s and @b Cs resides
/// @param [in, out] weights Pointer to the gsl_vector of the sample weights
/// @todo
/// Function name should be changed to "Shuffle"
void
altar::bayesian::COV::
rankAndShuffle(state_t & state, const cov::args & covargs, const gsl_permutation * perm) const
{
    // unpack the problem size
    const size_t samples = state.samples();
    const size_t parameters = state.parameters();

    // local variables
    size_t region, sample, samples_local, count;
    // coarse weights vector
    gsl_vector *w_coarse;
    // histogram
    gsl_histogram *h;

    // histogram ticks
    double *ticks; ///< @note <i>ticks</i> has to be double precision all the time, because gsl_histogram does not support single precision
    // counts
    //size_t * counts = new size_t[samples];
    gsl_vector * counts = gsl_vector_alloc(samples);

    // get the weight and weight_control
    const gsl_vector * w = covargs.w;
    const std::vector<cov::underflow_control_t> & wc = covargs.wc;
    // define the weight_hybrid gsl_vector
    gsl_vector * w_hybrid;
    // get the gsl_vector of coarse weights
    w_coarse = cov::coarseWeights(covargs);
    // get the iterator to the weight_control vector
    std::vector<cov::underflow_control_t>::const_iterator it;
    // get the number of weight regions
    const size_t regions = w_coarse->size;
    // sequentially loop over coarse_weight regions, except the last region
    size_t counts_left = samples;
    size_t sample_begin;
    for (region=0; region<regions-1; ++region)
    {
        it = covargs.wc.begin() + region;
        sample_begin = it->index_begin;
        samples_local = (it+1)->index_begin - sample_begin;
        // allocate w_hybrid
        w_hybrid = gsl_vector_alloc(samples_local+1); // hybrid w contains an extra element for the sum of the remaining weights
        // set the first samples_local elements of w_hybrid
        for (sample=0; sample<samples_local; ++sample) gsl_vector_set(w_hybrid, sample, gsl_vector_get(w, sample_begin+sample));
        // set the last element of w_hybrid
        gsl_vector_set(w_hybrid, samples_local, sumWeights(covargs, region+1));
        // shuffle inside this region
        ticks = new double[samples_local+2];
        ticks[0] = 0;
        for (sample=0; sample<samples_local+1; ++sample) ticks[sample+1] = ticks[sample] + gsl_vector_get(w_hybrid, sample);
        h = gsl_histogram_alloc(samples_local+1);
        gsl_histogram_set_ranges(h, ticks, samples_local+2);
        for (count=0; count<counts_left; ++count) gsl_histogram_increment(h, gsl_ran_flat(_rng, 0, ticks[samples_local+1]));
        for (sample=0; sample<samples_local; ++sample) gsl_vector_set(counts, sample_begin+sample, gsl_histogram_get(h, sample));
        counts_left = int(gsl_histogram_get(h, samples_local));
        // free locals
        delete [] ticks;
        gsl_histogram_free(h);
        gsl_vector_free(w_hybrid);
    } 
    // shuffle the last region
    samples_local = samples - wc.back().index_begin;
    w_hybrid = gsl_vector_alloc(samples_local);
    for (sample=0; sample<samples_local; ++sample) gsl_vector_set(w_hybrid, sample, gsl_vector_get(w, wc.back().index_begin+sample));
    ticks = new double[samples_local+1];
    ticks[0] = 0;
    for (sample=0; sample<samples_local; ++sample) ticks[sample+1] = ticks[sample] + gsl_vector_get(w_hybrid, sample);
    h = gsl_histogram_alloc(samples_local);
    gsl_histogram_set_ranges(h, ticks, samples_local+1);
    for (count=0; count<counts_left; ++count) gsl_histogram_increment(h, gsl_ran_flat(_rng, 0, ticks[samples_local]));
    for (sample=0; sample<samples_local; ++sample) gsl_vector_set(counts, wc.back().index_begin+sample, gsl_histogram_get(h, sample));
    // free locals
    delete [] ticks;
    gsl_histogram_free(h);
    gsl_vector_free(w_hybrid);
    

    // get the state vectors
    matrix_t * theta = state.theta();
    vector_t * prior = state.prior();
    vector_t * data = state.data();
    vector_t * posterior = state.posterior();

    // allocate duplicates
    matrix_t * thetaOld = gsl_matrix_TYPE_alloc(samples, parameters);
    vector_t * priorOld = gsl_vector_TYPE_alloc(samples);
    vector_t * dataOld = gsl_vector_TYPE_alloc(samples);
    vector_t * posteriorOld = gsl_vector_TYPE_alloc(samples);
    // make copies of all these
    gsl_matrix_TYPE_memcpy(thetaOld, theta);
    gsl_vector_TYPE_memcpy(priorOld, prior);
    gsl_vector_TYPE_memcpy(dataOld, data);
    gsl_vector_TYPE_memcpy(posteriorOld, posterior);

    // the number of samples we have processed
    size_t done = 0;
    // permute counts back to its corresponding order in state
    gsl_permutation * perm_inv = gsl_permutation_alloc(samples);
    gsl_permutation_inverse(perm_inv, perm);
    gsl_permute_vector(perm_inv, counts);
    // start shuffling the samples and likelihoods around
    for (size_t sample=0; sample<samples; ++sample) {
        // get its multiplicity
        count = static_cast<size_t>(gsl_vector_get(counts, sample));

        // get the row from the original samples set
        gsl_vector_TYPE_view row = gsl_matrix_TYPE_row(thetaOld, sample);
        // otherwise, duplicate this sample {count} number of times
        for (size_t dupl=0; dupl<count; ++dupl) {
            // by setting {count} consecutive rows of theta to this sample
            gsl_matrix_TYPE_set_row(theta, done, &row.vector);
            // and similarly for the log likelihoods
            gsl_vector_TYPE_set(prior, done, gsl_vector_TYPE_get(priorOld, sample));
            gsl_vector_TYPE_set(data, done, gsl_vector_TYPE_get(dataOld, sample));
            gsl_vector_TYPE_set(posterior, done, gsl_vector_TYPE_get(posteriorOld, sample));
            // update the {done} count
            done += 1;
        }
    }
    
    // free the temporaries
    gsl_vector_free(w_coarse);
    gsl_permutation_free(perm_inv);
    gsl_vector_free(counts);

    gsl_matrix_TYPE_free(thetaOld);
    gsl_vector_TYPE_free(priorOld);
    gsl_vector_TYPE_free(dataOld);
    gsl_vector_TYPE_free(posteriorOld);

    // all done
    return;
}

////////////////////////////////////////////////////////////////////////////////////////
// scale the weights and save the scaling information
// NOTE: this function assumes that llk in parameters has been sorted in descending order
// it also assumes that dbeta has been set in parameters
////////////////////////////////////////////////////////////////////////////////////////
void cov::scaleWeights(cov::args & p)
{
    const double dbeta = p.dbeta;

    // get the machine precision in log space (exponential entry)
    const double log_limit_min = log(std::numeric_limits<double>::min());
    const double log_limit_max = log(std::numeric_limits<double>::max());

    // local variables
    double w0_shift; // a positive value that shift w0 up to 1.
    double current_shift, new_shift; // a positive value that shift llk up
    double dbllk_shifted; // dbeta*shifted_llk
    underflow_control_t weight_control;

    // clear the wc first
    p.wc.clear();

    // save the scaling information of the first sample
    current_shift = w0_shift = std::fabs(dbeta*gsl_vector_get(p.llk,0));
    weight_control = {.index_begin=0, .scale=1.0}; // w0 serves as the standard for scaling
    p.wc.push_back(weight_control);

    // loop over each sample
    for (size_t i = 0; i < p.w->size; i++)
    {
        // get the shifted (dbeta*llk) value
        dbllk_shifted = dbeta*gsl_vector_get(p.llk,i) + current_shift;
        // if not shifted out of the underflow region...
        if (dbllk_shifted <= log_limit_min)
        {
            // set the new shift value that will shift the PREVIOUS sample (dbeta*llk) to be 0.
            new_shift = std::fabs(dbeta*gsl_vector_get(p.llk,i-1));
            // re-evaluate the shifted (dbeta*llk) value
            dbllk_shifted = dbeta*gsl_vector_get(p.llk,i) + new_shift;
            // if the shifted (dbeta*llk) is still not out of the underflow region
            if (dbllk_shifted <= log_limit_min)
            {
                ///< @note should we proceed?
                // send warning/error if still not shifted out of the underflow region
                // altar::utils::Warning("Sample cannot to shifted out of the underflow region!");
                // set the new shift value that will shift the CURRENT sample (dbeta*llk) to be 0.
                new_shift = std::fabs(dbeta*gsl_vector_get(p.llk,i));
                // and re-evaluate the shifted (dbeta*llk) value
                dbllk_shifted = dbeta*gsl_vector_get(p.llk,i) + new_shift;
                // and save the scaling information
                weight_control = {.index_begin=i, .scale=std::exp(-(new_shift-current_shift))};
                p.wc.push_back(weight_control);
            }
            // otherwise (if the shifted (dbeta*llk) is successfully out of the underflow region)
            else
            {
                // if there was only 1 sample processed, remove the previous weight_control 
                if (i==1) p.wc.pop_back();
                // reset the previous sample weight to be 1.0
                gsl_vector_set(p.w, i-1, 1.0);
                // and save the scaling information
                weight_control = {.index_begin=i-1, .scale=std::exp(-(new_shift-current_shift))};
                p.wc.push_back(weight_control);
            }
            // save the new shift value
            current_shift = new_shift;
        } 
        // set the scaled weight
        gsl_vector_set(p.w, i, std::exp(dbllk_shifted));
    }
}

////////////////////////////////////////////////////////////////////////////////////////
// return a gsl_vector that contains the fine weight information for each sample
////////////////////////////////////////////////////////////////////////////////////////
gsl_vector * cov::fineWeights(const args & p)
{
    // get the system size
    const size_t samples = p.w->size;

    // allocate the fine_weight vector
    gsl_vector * fine_weight= gsl_vector_alloc(samples);

    // local variables
    size_t sample, samples_local, index_begin;

    // do the work
    double scale=1.0;
    for (std::vector<underflow_control_t>::const_iterator it=p.wc.begin(); it!=p.wc.end(); ++it)
    {
        index_begin = it->index_begin;
        if ((it+1)!=p.wc.end()) samples_local = (it+1)->index_begin - index_begin;
        else samples_local = samples - index_begin;
        scale *= it->scale;
        for (sample=0; sample<samples_local; ++sample)
        {
            gsl_vector_set(fine_weight, index_begin+sample, gsl_vector_get(p.w, index_begin+sample)*scale);
        } 
    } 

    // return the fine_weight gsl_vector
    return fine_weight;
}

////////////////////////////////////////////////////////////////////////////////////////
// return a gsl_vector that contains the unscaled coarse (sum) weight information for each sub_region
////////////////////////////////////////////////////////////////////////////////////////
gsl_vector * cov::coarseWeights(const args & p)
{
    // get the system size
    const size_t samples = p.w->size;

    // allocate the coarse_weight vector
    const size_t region = p.wc.size();
    gsl_vector * coarse_weight = gsl_vector_alloc(region);

    // assign sum_of_weights for each region
    size_t count=0;
    size_t index_begin;
    gsl_vector_view w_sub;
    for (std::vector<underflow_control_t>::const_iterator it=p.wc.begin(); it!=p.wc.end(); ++it)
    {
        index_begin = it->index_begin;
        if ((it+1)!=p.wc.end()) w_sub = gsl_vector_subvector(p.w, index_begin, ((it+1)->index_begin - index_begin) );
        else w_sub = gsl_vector_subvector(p.w, index_begin, (samples - index_begin) );
        gsl_vector_set(coarse_weight, count, gsl_blas_dasum(&w_sub.vector));
        ++ count;
    } 

    // return the coarse_weight gsl_vector
    return coarse_weight;
}


////////////////////////////////////////////////////////////////////////////////////////
// get the sum of the scaled-weights beyond a given region (including this region)
////////////////////////////////////////////////////////////////////////////////////////
double cov::sumWeights(const args & p, const int region_start)
{
    // get the system size
    const size_t samples = p.w->size;
    // and number of weight regions
    const size_t regions = p.wc.size();

    // get the coarse weights
    gsl_vector * w_coarse = cov::coarseWeights(p);

    // local variables
    double coarse_weight;

    // assign sum_of_weights for each region
    size_t count=0;
    double sum=0.0;
    for (int i=regions-1; i>=region_start; --i)
    {
        coarse_weight = gsl_vector_get(w_coarse, i);
        sum += coarse_weight;
        sum *= p.wc.at(i).scale;
    } 

    // free locals
    gsl_vector_free(w_coarse);

    // return the sum of the scaled weights
    return sum;
}

////////////////////////////////////////////////////////////////////////////////////////
// get the sum of the squared scaled-weights
// it is only being used by cov::cov function, and it does not have the same prototype as sumWeights
////////////////////////////////////////////////////////////////////////////////////////
double cov::sumWeightsSquared(const args & p)
{
    // get the system size
    const size_t samples = p.w->size;
    // and number of weight regions
    const size_t regions = p.wc.size();

    // assign sum_of_weights for each region
    double scale;
    double sum_squared=0.0;
    int index_begin;
    gsl_vector_view w_sub;
    std::vector<underflow_control_t>::const_iterator it;
    for (int i=regions-1; i>=0; --i)
    {
        it = p.wc.begin()+i;
        index_begin = it->index_begin;
        if ((it+1)!=p.wc.end()) w_sub = gsl_vector_subvector(p.w, index_begin, ((it+1)->index_begin - index_begin) );
        else w_sub = gsl_vector_subvector(p.w, index_begin, (samples - index_begin) );
        sum_squared += pow((gsl_blas_dnrm2(&w_sub.vector)),2.);
        scale = p.wc.at(i).scale;
        sum_squared *= (scale*scale);
    }

    // return the sum of the squred scaled-weights
    return sum_squared;
}


////////////////////////////////////////////////////////////////////////////////////////
// calculate the cov value
// NOTE: this function assumes that dataLLK has been sorted in descending order
////////////////////////////////////////////////////////////////////////////////////////
double cov::cov(double dbeta, void * parameters)
{
    // get my auxiliary parameters
    args & p = *static_cast<cov::args *>(parameters);
    // store dbeta
    p.dbeta = dbeta;

    /*
    /////////////////////////////////
    // initialize {w}
    gsl_vector * w = p.w;
    double median = gsl_stats_max(w->data,w->stride,w->size);
    for (size_t i = 0; i < p.w->size; i++) {
        gsl_vector_set(p.w, i, std::exp(dbeta * (gsl_vector_get(p.llk, i)  - median)));
    }

    // normalize
    double wsum = p.w->size * gsl_stats_mean(p.w->data, p.w->stride, p.w->size);
    gsl_vector_scale(p.w, 1/wsum);

    // compute the COV
    double mean = gsl_stats_mean(p.w->data, p.w->stride, p.w->size);
    double sdev = gsl_stats_sd(p.w->data, p.w->stride, p.w->size);
    /////////////////////////////////
    */

    // std::cout << " ** altar.beta: dbeta = " << dbeta << std::endl;
    
    // initialize and scale {w}
    scaleWeights(p);

    // compute the COV
    double cov;
    double S=0., SS=0.;
    const size_t N = p.w->size;
    S = sumWeights(p);
    SS = sumWeightsSquared(p);
    cov = sqrt((N*N*SS-N*S*S)/(N-1)) / S;

    // store it
    p.cov = cov;

    // if the COV is not well-defined
    if (gsl_isinf(cov) || gsl_isnan(cov)) {
        // set our metric to some big value
        p.metric = 1e100;
    } else {
        // otherwise
        p.metric = gsl_pow_2(cov - p.target);
    }
    
    // and return
    return p.metric;
}

// end of file
