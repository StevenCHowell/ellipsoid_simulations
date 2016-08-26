// -*- C++ -*-
// 
// Hailiang Zhang
// california institute of technology
// 


// for the build system
#include <portinfo>

// get my declaration
#include "cudaMetropolis.h"

// my dependencies
#include <altar/problem/cudaProblem.h>
#include <altar/problem/cudaCoolingStep.h>


/// @par Main functionality
/// Random walk on GPU
/// @param [in, out] state The CPU state where @b M, <i>LLK</i>'s and @b Cs resides
/// @param [in, out] cuda_state The GPU state where @b M, <i>LLK</i>'s and @b Cs resides
/// @param [in] problem The problem entity on GPU (including priors, models, data)
altar::bayesian::cuda::cudaMetropolis::stats_t
altar::bayesian::cuda::cudaMetropolis::
sample(state_t & state, cuda_state_t & cuda_state, const cuda_problem_t & problem) const
{
    // some profiling
    bool P_details=false; 
    struct timeval t1, t2, t1all, t2all;
    bool P_altar_profiling, P_altar_metropolis_detail_profiling;
    if (getenv("ALTAR_PROFILING") && strcmp(getenv("ALTAR_PROFILING"),"true")==0) P_altar_profiling = true;
    if (getenv("ALTAR_METROPOLIS_DETAIL_PROFILING") && strcmp(getenv("ALTAR_METROPOLIS_DETAIL_PROFILING"),"true")==0) P_altar_metropolis_detail_profiling = true;

    // get the problem sizes
    const size_t Ns = state.samples();
    const size_t NMparam = state.parameters();


    // get theta
    matrix_t * const thetaT = gsl_matrix_TYPE_alloc(NMparam, Ns);
    // get data covariance matrix
    matrix_t * const sigma = state.sigma();

    // get the current gpu number
    const int ngpu = cuda_state.ngpu();
    const int rank= cuda_state.rank();

    // get sigma_chol on GPU
    TYPE * const gCmi = cuda_state.gpu_sigma();
    // get M on GPU
    TYPE * const gM = cuda_state.gpu_M();
    // get M_candidate on GPU
    TYPE * const gM_candidate = cuda_state.gpu_M_candidate();
    // get LLK on GPU
    TYPE * const gLLKprior = cuda_state.gpu_prior();
    TYPE * const gLLKdata = cuda_state.gpu_data();
    TYPE * const gLLKpost = cuda_state.gpu_posterior();
    // get candidate LLK on GPU
    TYPE * const gLLKprior_candidate = cuda_state.gpu_prior_candidate();
    TYPE * const gLLKdata_candidate = cuda_state.gpu_data_candidate();
    TYPE * const gLLKpost_candidate = cuda_state.gpu_posterior_candidate();
    // get pnacc on GPU
    int * const gpnacc = cuda_state.gpu_pnacc();
    // get random number array for walk on GPU
    TYPE * const gR_walk = cuda_state.gpu_R_walk();
    // get random number array for Metropolis selection on GPU
    TYPE * const gR_select = cuda_state.gpu_R_select();

    // get the rejection flags on GPU
    int * gpu_rejection_flags = cuda_state.gpu_rejection_flags();

    // allocate room for my Cholesky decomposed covariance matrix
    matrix_t * sigma_chol = gsl_matrix_TYPE_alloc(NMparam, NMparam);
    // make a copy
    gsl_matrix_TYPE_memcpy(sigma_chol, sigma);
    //altar::problem::util::mprintf(sigma_chol, "Sigma original:", 15, 15);

    // scale it
    gsl_matrix_TYPE_scale(sigma_chol, _scaling);
    //altar::problem::util::mprintf(sigma_chol, "Sigma before conditioning:", 15, 15);

    // Cholesky decompose
#ifdef USE_DOUBLE
    gsl_linalg_cholesky_decomp(sigma_chol);
#else
    gsl_matrix * sigma_chol_double = altar::utils::gsl_matrix_double_from_float(sigma_chol);
    gsl_linalg_cholesky_decomp(sigma_chol_double);
    altar::utils::gsl_matrix_set_float_from_double(sigma_chol, sigma_chol_double);
    gsl_matrix_free(sigma_chol_double);
#endif
    //altar::problem::util::mprintf(sigma_chol, "Sigma after decomposition:", 15, 15);

    // copy back
    gsl_matrix_TYPE_memcpy(sigma, sigma_chol);
    //altar::problem::util::mprintf(state.theta(), "theta:", 15, 15);

    // write to externel files
    static int count=0;
    size_t i, j;
    std::string cuda_statistics_dir = "cuda_statistics/"; // it's important that directory name is postfixed with a slash
    if (!opendir(cuda_statistics_dir.c_str())) mkdir(cuda_statistics_dir.c_str(), S_IRWXU);
    //else system((std::string("exec rm -rf ")+cuda_statistics_dir+std::string("/*")).c_str());
    std::stringstream ss_overall, ss_details;
    ss_overall<<cuda_statistics_dir<<"overall_rank_"<<rank<<"_gpu_"<<ngpu<<".txt";
    ss_details<<cuda_statistics_dir<<"details_rank_"<<rank<<"_gpu_"<<ngpu<<".txt";
    std::ofstream fout_overall, fout_details;
    if (!count)
    {
        fout_overall.open(ss_overall.str().c_str());
        fout_details.open(ss_details.str().c_str());
    }
    else
    {
        fout_overall.open(ss_overall.str().c_str(), std::ios::app);
        fout_details.open(ss_details.str().c_str(), std::ios::app);
    }
    fout_overall<<std::endl<<"================================"<<std::endl<<"Iteration: "<<count<<std::endl;
    fout_details<<std::endl<<"================================"<<std::endl<<"Iteration: "<<count<<std::endl;
    fout_overall<<"Beta: "<<state.beta()<<std::endl;
    fout_overall.setf(std::ios::fixed);
    fout_overall.precision(3);
    fout_overall<<"Starting theta:"<<std::endl;
    for (i=0; i<(Ns<10?Ns:10); i++) { for(j=0; j<(NMparam<18?NMparam:18); j++) fout_overall<<std::setw(8)<<gsl_matrix_TYPE_get(state.theta(),i,j)<<" "; fout_overall<<std::endl;}
    fout_overall<<"Starting LLKs:"<<std::endl;
    fout_overall<<gsl_stats_TYPE_mean(state.prior()->data, 1, Ns)<<std::endl;
    for (i=0; i<(Ns<10?Ns:10); i++) { fout_overall<<std::setw(12)<<gsl_vector_TYPE_get(state.prior(),i)<<" "; } fout_overall<<"..."<<std::endl;
    fout_overall<<"..."; for (i=0; i<(Ns<10?Ns:10); i++) { fout_overall<<std::setw(12)<<gsl_vector_TYPE_get(state.prior(),Ns-1-i)<<" "; } fout_overall<<std::endl;
    fout_overall<<gsl_stats_TYPE_mean(state.data()->data, 1, Ns)<<std::endl;
    for (i=0; i<(Ns<10?Ns:10); i++) { fout_overall<<std::setw(12)<<gsl_vector_TYPE_get(state.data(),i)<<" "; } fout_overall<<"..."<<std::endl;
    fout_overall<<"..."; for (i=0; i<(Ns<10?Ns:10); i++) { fout_overall<<std::setw(12)<<gsl_vector_TYPE_get(state.data(),Ns-1-i)<<" "; } fout_overall<<std::endl;
    fout_overall<<gsl_stats_TYPE_mean(state.posterior()->data, 1, Ns)<<std::endl;
    for (i=0; i<(Ns<10?Ns:10); i++) { fout_overall<<std::setw(12)<<gsl_vector_TYPE_get(state.posterior(),i)<<" "; } fout_overall<<"..."<<std::endl;
    fout_overall<<"..."; for (i=0; i<(Ns<10?Ns:10); i++) { fout_overall<<std::setw(12)<<gsl_vector_TYPE_get(state.posterior(),Ns-1-i)<<" "; } fout_overall<<std::endl<<std::endl;


    ///////////////////////////////////////////////////
    // distribute theta, Cmi, and llk to gpu    
    ///////////////////////////////////////////////////
    cuda_state.distributeGPU(state);

    //altar::utils::cuda::print_gpu_matrix(fout_details, gCmi, NMparam, NMparam, 10, 18, "gCmi:");
    //altar::utils::cuda::print_gpu_matrix(fout_details, gM, Ns, NMparam, 10, 18, "gM:");//,"%8.3f");

    // set gpu_nacc to be 0
    CALL_CUDA(cudaMemset(gpnacc,0,sizeof(int)));

    // initialize the counts
    size_t acceptedSamples = 0;
    size_t rejectedSamples = 0;
    

    // Hack the starting time
    if (P_altar_profiling) gettimeofday(&t1all,NULL);
    if (P_altar_metropolis_detail_profiling) gettimeofday(&t1,NULL);

    ///////////////////////////////////////////////////
    // walk the chains
    ///////////////////////////////////////////////////
    for (size_t link=0; link<_steps; ++link)
    {

        // random displacement
#ifdef USE_DOUBLE 
        CALL_CURAND(curandGenerateNormalDouble(_curand_gen, gR_walk, Ns*NMparam, 0, 1));
#else
        CALL_CURAND(curandGenerateNormal(_curand_gen, gR_walk, Ns*NMparam, 0, 1));
#endif

        /*
        // ... this is for debugging purpose only
        matrix_t * delta = gsl_matrix_TYPE_alloc(Ns, NMparam);
        // fill it
        for (size_t i=0; i<Ns; ++i) {
            for (size_t j=0; j<NMparam; ++j) {
                gsl_matrix_TYPE_set(delta, i, j, gsl_ran_ugaussian(_rng));
            }   
        }   
        matrix_t * deltaT = gsl_matrix_TYPE_alloc(NMparam, Ns);
        gsl_matrix_TYPE_transpose_memcpy(deltaT, delta);
        CALL_CUDA(cudaMemcpy(gR_walk, deltaT->data, Ns*NMparam*sizeof(TYPE), cudaMemcpyHostToDevice));
        */

        if (P_details) altar::utils::cuda::print_gpu_matrix(fout_details, gR_walk, Ns, NMparam, 10, 18, "gR_walk:");

#ifdef USE_DOUBLE 
        CALL_CUBLAS(cublasDtrmm(_cublas_handle, CUBLAS_SIDE_RIGHT, CUBLAS_FILL_MODE_LOWER, CUBLAS_OP_T, CUBLAS_DIAG_NON_UNIT, Ns, NMparam, &const_1, gCmi, NMparam, gR_walk, Ns, gM_candidate, Ns));//May have problem
        CALL_CUBLAS(cublasDaxpy(_cublas_handle, NMparam*Ns, &const_1, gM, 1, gM_candidate, 1));        
#else
        CALL_CUBLAS(cublasStrmm(_cublas_handle, CUBLAS_SIDE_RIGHT, CUBLAS_FILL_MODE_LOWER, CUBLAS_OP_T, CUBLAS_DIAG_NON_UNIT, Ns, NMparam, &const_1, gCmi, NMparam, gR_walk, Ns, gM_candidate, Ns));//May have problem
        CALL_CUBLAS(cublasSaxpy(_cublas_handle, NMparam*Ns, &const_1, gM, 1, gM_candidate, 1));        
#endif
        if (P_details) altar::utils::cuda::print_gpu_matrix(fout_details, gM, Ns, NMparam, 10, 18, "gM:");//,"%8.3f");
        if (P_details) altar::utils::cuda::print_gpu_matrix(fout_details, gM_candidate, Ns, NMparam, 10, 18, "gM_candidate:");//,"%8.3f");
        if (P_details) altar::utils::cuda::print_gpu_matrix(fout_details, cuda_state.gpu_M_candidate_queued(), Ns, NMparam, 10, 18, "gM_candidate_queued:");//,"%8.3f");
        if (P_altar_metropolis_detail_profiling) altar::utils::print_profiling(rank,t1, t2, (char*)"Metropolis-displacement");

        // ask the problem to verify it's bad move or not
        problem.verify(cuda_state, gpu_rejection_flags);
        if (P_details) altar::utils::cuda::print_gpu_matrix(fout_details, gpu_rejection_flags, Ns, 1, 50, 1, "rejection flags:", 2, 0);
        if (P_altar_metropolis_detail_profiling) altar::utils::print_profiling(rank,t1, t2, (char*)"Metropolis-verify");

        // ask the problem to compute the log likelihoods of the candidate
        problem.likelihoods(cuda_state);
        if (P_details) altar::utils::cuda::print_gpu_vector(fout_details, gLLKprior_candidate, Ns, 10, "gLLKprior_candidate:");
        if (P_details) altar::utils::cuda::print_gpu_vector(fout_details, gLLKdata_candidate, Ns, 10, "gLLKdata_candidate:");
        if (P_details) altar::utils::cuda::print_gpu_vector(fout_details, gLLKpost_candidate, Ns, 10, "gLLKpost_candidate:");
        if (P_altar_metropolis_detail_profiling) altar::utils::print_profiling(rank,t1, t2, (char*)"Metropolis-likelihoods");


        // Metropolis screening
#ifdef USE_DOUBLE 
        CALL_CURAND(curandGenerateUniformDouble(_curand_gen, gR_select, Ns));
#else
        CALL_CURAND(curandGenerateUniform(_curand_gen, gR_select, Ns));
#endif
        if (P_details) altar::utils::cuda::print_gpu_vector(fout_details, gR_select, Ns, 10, "gR_select:");
        wrapperCudaUpdateSample(Ns, NMparam, gpu_rejection_flags, gR_select, gpnacc,
                                gM, gLLKprior, gLLKdata, gLLKpost,
                                gM_candidate, gLLKprior_candidate, gLLKdata_candidate, gLLKpost_candidate);
        if (P_details) altar::utils::cuda::print_gpu_matrix(fout_details, gM, Ns, NMparam, 10, 18, "gM after update:");//,"%8.3f");
        if (P_details) altar::utils::cuda::print_gpu_vector(fout_details, gpnacc, 1, 1, "gpnacc:");
        if (P_altar_metropolis_detail_profiling) altar::utils::print_profiling(rank,t1, t2, (char*)"Metropolis-update");
    }

    if (P_altar_profiling) altar::utils::print_profiling(rank,t1all, t2all, (char*)"Metropolis");

    ///////////////////////////////////////
    // GPU collection
    ///////////////////////////////////////
    cuda_state.collectGPU(state);

    // free the temporaries
    gsl_matrix_TYPE_free(sigma_chol);
    gsl_matrix_TYPE_free(thetaT);


    // get the number of the accepted samples
    int * pacceptedSamples;
    CALL_CUDA(cudaMallocHost((void**)&pacceptedSamples, sizeof(int)));
    CALL_CUDA(cudaMemcpy(pacceptedSamples, gpnacc, sizeof(int), cudaMemcpyDeviceToHost));
    acceptedSamples = *pacceptedSamples;
    CALL_CUDA(cudaFreeHost(pacceptedSamples));
    rejectedSamples = _steps*Ns- acceptedSamples;
    // printf("GPU #%d Accepted: %d, Rejected: %d, Acceptance Rata: %f\n", ngpu, acceptedSamples, rejectedSamples, acceptedSamples/TYPE(acceptedSamples+rejectedSamples));

    /*
    // write to externel files
    fout_overall<<"\nAccepted: "<<acceptedSamples<<" Rejected: "<<rejectedSamples<<" Acceptance rate: "<<acceptedSamples/TYPE(acceptedSamples+rejectedSamples)<<std::endl;
    fout_overall<<"\nEnding theta:"<<std::endl;
    for (i=0; i<(Ns<10?Ns:10); i++) { for(j=0; j<(NMparam<18?NMparam:18); j++) fout_overall<<std::setw(8)<<gsl_matrix_TYPE_get(state.theta(),i,j)<<" "; fout_overall<<std::endl;}
    fout_overall<<"Ending LLKs:"<<std::endl;
    fout_overall<<gsl_stats_TYPE_mean(state.prior()->data, 1, Ns)<<std::endl;
    for (i=0; i<(Ns<10?Ns:10); i++) { fout_overall<<std::setw(12)<<gsl_vector_TYPE_get(state.prior(),i)<<" "; } fout_overall<<"..."<<std::endl;
    fout_overall<<"..."; for (i=0; i<(Ns<10?Ns:10); i++) { fout_overall<<std::setw(12)<<gsl_vector_TYPE_get(state.prior(),Ns-1-i)<<" "; } fout_overall<<std::endl;
    fout_overall<<gsl_stats_TYPE_mean(state.data()->data, 1, Ns)<<std::endl;
    for (i=0; i<(Ns<10?Ns:10); i++) { fout_overall<<std::setw(12)<<gsl_vector_TYPE_get(state.data(),i)<<" "; } fout_overall<<"..."<<std::endl;
    fout_overall<<"..."; for (i=0; i<(Ns<10?Ns:10); i++) { fout_overall<<std::setw(12)<<gsl_vector_TYPE_get(state.data(),Ns-1-i)<<" "; } fout_overall<<std::endl;
    fout_overall<<gsl_stats_TYPE_mean(state.posterior()->data, 1, Ns)<<std::endl;
    for (i=0; i<(Ns<10?Ns:10); i++) { fout_overall<<std::setw(12)<<gsl_vector_TYPE_get(state.posterior(),i)<<" "; } fout_overall<<"..."<<std::endl;
    fout_overall<<"..."; for (i=0; i<(Ns<10?Ns:10); i++) { fout_overall<<std::setw(12)<<gsl_vector_TYPE_get(state.posterior(),Ns-1-i)<<" "; } fout_overall<<std::endl<<std::endl;
*/

    // close externel files
    fout_overall.close();
    fout_details.close();
    ++count;


    // and return
    return stats_t(acceptedSamples, rejectedSamples);
}


/// @par Main functionality
/// destructor
/// @note
/// The cublas handle and cuda random generator are destroyed here
altar::bayesian::cuda::cudaMetropolis::
~cudaMetropolis() 
{
    CALL_CUBLAS(cublasDestroy(_cublas_handle));
    CALL_CURAND(curandDestroyGenerator(_curand_gen));
}

// implementation details


// end of file
