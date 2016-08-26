// -*- C++ -*-
// 
// michael a.g. aïvázis
// california institute of technology
// (c) 2010-2013 all rights reserved
// 


// for the build system
#include <portinfo>

// get my declarations
#include "AnnealingMethod.h"
// and my dependencies
#include "COV.h"

// externels
#include "H5Cpp.h"


/// @par Main functionality
/// Update model data (eg gf, d) for Cp implementation
/// @param [in, out] annealer Host of the problem entity where priors, models and data reside
/// @return True if Cp implemented, false if Cp is not implemented
bool
altar::bayesian::AnnealingMethod::
updateModelData(annealer_t &)
{
    return false;
}


/// @par Main functionality
/// Start the annealing process from scratch<br>
/// @note
/// This is just a dummy method to set the iteration number to be 0<br>
/// The actually implementation is done by its derived classes
altar::bayesian::AnnealingMethod &
altar::bayesian::AnnealingMethod::
start(annealer_t &, const std::string &)
{
    // reset my iteration count
    _iteration = 0;
    // all done
    return *this;
}


/// @par Main functionality
/// Start the annealing process from a previous run<br>
/// @param [in, out] annealer Host of the problem entity where priors, models and data reside
/// @param [in] iteration The iteration number to be restarted
altar::bayesian::AnnealingMethod &
altar::bayesian::AnnealingMethod::
restart(annealer_t & annealer, const size_t iteration)
{
    // reset my iteration count
    _iteration = iteration;

    // make sure the output directory "results" exist
    std::string results_dir = "results/"; // it's important that directory name is postfixed with a slash
    if (!opendir(results_dir.c_str())) std::cerr<<"Folder: "<<results_dir<<" does not exit!"<<std::endl;

    // get beta and scaling from BetaStatistics.txt
    TYPE beta, scaling;
    std::string file_statistics = results_dir+"BetaStatistics.txt";
    _readBetastatistics(file_statistics, iteration, beta, scaling);

    // set my beta
    scheduler_t & scheduler = annealer.scheduler();
    scheduler.beta(beta);
    _state.beta(beta);
    
    // set my scaling factor
    sampler_t & sampler = annealer.sampler();
    sampler._scaling = scaling; // Hailiang made them friends

    // all done
    return *this;
}

/// @par Main functionality
/// Notify me when annealing is finished
/// @param [in] annealer Host of the problem entity where priors, models and data reside
/// @note
/// This is just a dummy method to set the iteration number to be 0<br>
/// The actually implementation is done by its derived classes
altar::bayesian::AnnealingMethod &
altar::bayesian::AnnealingMethod::
finish(annealer_t &)
{
    // all done
    return *this;
}

/// @par Main functionality
/// Push my state forward along the cooling schedule
/// @param [in, out] annealer Host of the problem entity where priors, models and data reside
altar::bayesian::AnnealingMethod &
altar::bayesian::AnnealingMethod::
cool(annealer_t & annealer)
{
    // get the scheduler
    annealer_t::scheduler_t & scheduler = annealer.scheduler();
    // ask it to update my state
    scheduler.update(_state);
    // update my iteration counter
    iterate();
    // all done
    return *this;
}

/// @par Main functionality
/// Increment the iteration number by 1
void
altar::bayesian::AnnealingMethod::
iterate()
{
    _iteration += 1;
}

/// @par Main functionality
/// Re-sample the posterior distribution
/// @param [in, out] annealer Host of the problem entity where priors, models and data reside
altar::bayesian::AnnealingMethod::stats_t
altar::bayesian::AnnealingMethod::
resample(annealer_t & annealer)
{
    // get the problem
    annealer_t::problem_t & problem = annealer.problem();
    // and the sampler
    annealer_t::sampler_t & sampler = annealer.sampler();

    // ask it to sampler the posterior pdf and return the sampling statistics
    return sampler.sample(_state, problem);
}


/// @par Main functionality
/// Analyze the acceptance statistics and take the problem state to the end of the current annealing step
/// @param [in, out] annealer Host of the problem entity where priors, models and data reside
/// @param [in, out] stats Alias for (acceptedSamples, rejectedSamples)
altar::bayesian::AnnealingMethod &
altar::bayesian::AnnealingMethod::
equilibrate(annealer_t & annealer, const annealer_t::stats_t & stats)
{
    // get the sampler
    annealer_t::sampler_t & sampler = annealer.sampler();
    // ask it to equilibrate
    sampler.equilibrate(stats);

    // all done
    return *this;
}

/// @par Main functionality
/// Status report
/// @param [in, out] annealer Host of the problem entity where priors, models and data reside
altar::bayesian::AnnealingMethod &
altar::bayesian::AnnealingMethod::
report(annealer_t & annealer)
{
    // if my info channel is active
    if (_info) {
        // ask for a report
        _report(annealer);
        printf("Theta_mean: \n");
        gsl_vector_TYPE * v = gsl_vector_TYPE_calloc(_state.samples());
        for (size_t i=0; i<_state.parameters(); i++)
        {
            gsl_matrix_TYPE_get_col(v,_state.theta(),i);
            printf("%f ",gsl_stats_TYPE_mean(v->data,v->stride,v->size));
        }
        printf("\n");
        printf("Theta_sd: \n");
        for (size_t i=0; i<_state.parameters(); i++)
        {
            gsl_matrix_TYPE_get_col(v,_state.theta(),i);
            printf("%f ",gsl_stats_TYPE_sd(v->data,v->stride,v->size));
        }
        printf("\n");
        //printf("\nInitial LLK_prior: \n"); for (size_t i=0; i<_state.samples(); i++) printf("%f ",gsl_vector_TYPE_get(_state.prior(),i));
        //printf("\nInitial LLK_data: \n"); for (size_t i=0; i<_state.samples(); i++) printf("%f ",gsl_vector_TYPE_get(_state.data(),i));
        //printf("\nInitial LLK_posterior: \n"); for (size_t i=0; i<_state.samples(); i++) printf("%f ",gsl_vector_TYPE_get(_state.posterior(),i));
        // flush
        _info << pyre::journal::endl;
    }

    // all done
    return *this;
}


/// @par Main functionality
/// Write the annealing progress to the "results" folder
/// @param [in] annealer Host of the problem entity where priors, models and data reside
/// @param [in, out] stats Alias for (acceptedSamples, rejectedSamples)
/// @param [in] restart The restart iteration number
altar::bayesian::AnnealingMethod &
altar::bayesian::AnnealingMethod::
write_out(annealer_t & annealer, const stats_t & stats, const bool flag_detail, const size_t restart)
{
    // make sure the output directory "results" exist
    std::string results_dir = "results/"; // it's important that directory name is postfixed with a slash
    if (!opendir(results_dir.c_str())) mkdir(results_dir.c_str(), S_IRWXU);
    //else system((std::string("exec rm -rf ")+results_dir).c_str());


    // write out staticstics
    std::string file_statistics = results_dir+"BetaStatistics.txt";
    std::fstream fp_statistics;
    if (!_iteration) //if calling it for the first time
    {
        fp_statistics.open(file_statistics.c_str(), std::ios::out);
        fp_statistics<<"#Step  |  Scaling  |  COV  |  Beta  |  Accepted  |  Rejected \n";
    }
    else
    {
        if (restart && restart==_iteration)
        {
            fp_statistics.open(file_statistics.c_str(), std::ios::in);
            std::stringstream ss;
            std::string line;
            std::vector<std::string> lines;
            size_t iteration;
            while (std::getline(fp_statistics, line))
            {
                ss.str(line);
                if (ss>>iteration, iteration>restart) break; //move the file pointer to the iteration after restart number
                lines.push_back(line);
                ss.clear();
            }
            ss.clear();
            fp_statistics.close();
            fp_statistics.open(file_statistics.c_str(), std::ios::out);
            for (std::vector<std::string>::iterator it = lines.begin(); it!=lines.end(); ++it) fp_statistics<<*it<<std::endl;
            fp_statistics.close();
        }
        else fp_statistics.open(file_statistics.c_str(), std::ios::out|std::ios::app);
    }
        
    fp_statistics<<_iteration<<" "<<annealer.sampler()._scaling<<" "<<annealer.scheduler().cov()<<" "<<annealer.scheduler().beta()<<" "
                 <<std::get<0>(stats)<<" "<<std::get<1>(stats)<<std::endl;
    fp_statistics.close();
    
    if (flag_detail)
    {
        // set the hdf5 file name
        std::stringstream ss;
        ss<<results_dir<<"step_"<<std::setw(3)<<std::setfill('0')<<_iteration<<".h5";

        // create the hdf5 file
        H5::H5File * h5file = new H5::H5File(ss.str(),H5F_ACC_TRUNC);

        // hdf5 stuff
        hsize_t * fdim = new hsize_t(2);
        H5::DataSpace fspace;
        H5::DataSet dataset;
        H5::Attribute attribute;
        H5::DataSpace attribute_space = H5::DataSpace(H5S_SCALAR);
        H5::StrType attribute_strtype = H5::StrType(H5::PredType::C_S1, 256);

        // write out COV 
        fdim[0]=_state.parameters();
        fdim[1]=_state.parameters();
        fspace = H5::DataSpace(2, fdim);
        dataset = H5::DataSet(h5file->createDataSet("Covariance", H5::PredType::NATIVE_TYPE, fspace));
        dataset.write(state().sigma_const()->data, H5::PredType::NATIVE_TYPE, fspace);
        attribute = dataset.createAttribute(H5std_string("Help"), attribute_strtype, attribute_space);
        attribute.write(attribute_strtype, H5std_string("Covariance matrix of the current sample set (parameters * parameters)"));

        // write out theta
        fdim[0]=_state.samples();
        fdim[1]=_state.parameters();
        fspace = H5::DataSpace(2, fdim);
        dataset = H5::DataSet(h5file->createDataSet("Sample Set", H5::PredType::NATIVE_TYPE, fspace));
        dataset.write(state().theta_const()->data, H5::PredType::NATIVE_TYPE, fspace);
        attribute = dataset.createAttribute(H5std_string("Help"), attribute_strtype, attribute_space);
        attribute.write(attribute_strtype, H5std_string("Set of samples of model parameters (samples * parameters)"));

        // write out prior llk
        fdim[0]=_state.samples();
        fspace = H5::DataSpace(1, fdim);
        dataset = H5::DataSet(h5file->createDataSet("Prior Log-likelihood", H5::PredType::NATIVE_TYPE, fspace));
        dataset.write(state().prior_const()->data, H5::PredType::NATIVE_TYPE, fspace);
        attribute = dataset.createAttribute(H5std_string("Help"), attribute_strtype, attribute_space);
        attribute.write(attribute_strtype, H5std_string("Natural log of prior likelihoods (samples)"));


        // write out data llk
        fdim[0]=_state.samples();
        fspace = H5::DataSpace(1, fdim);
        dataset = H5::DataSet(h5file->createDataSet("Data Log-likelihood", H5::PredType::NATIVE_TYPE, fspace));
        dataset.write(state().data_const()->data, H5::PredType::NATIVE_TYPE, fspace);
        attribute = dataset.createAttribute(H5std_string("Help"), attribute_strtype, attribute_space);
        attribute.write(attribute_strtype, H5std_string("Natural log of data likelihoods (samples)"));

        // write out posterior llk
        fdim[0]=_state.samples();
        fspace = H5::DataSpace(1, fdim);
        dataset = H5::DataSet(h5file->createDataSet("Posterior Log-likelihood", H5::PredType::NATIVE_TYPE, fspace));
        dataset.write(state().posterior_const()->data, H5::PredType::NATIVE_TYPE, fspace);
        attribute = dataset.createAttribute(H5std_string("Help"), attribute_strtype, attribute_space);
        attribute.write(attribute_strtype, H5std_string("Natural log of posterior likelihoods (samples)"));

        // clean hdf5 stuff
        delete h5file;
        delete fdim;
    }

    // all done
    return *this;
}


/// @par Main functionality
/// destructor
altar::bayesian::AnnealingMethod::
~AnnealingMethod() 
{}

/// @par Main functionality
/// status report engine
/// @param [in, out] annealer Host of the problem entity where priors, models and data reside
altar::bayesian::AnnealingMethod &
altar::bayesian::AnnealingMethod::
_report(annealer_t & annealer)
{
    _info
        << pyre::journal::at(__HERE__)
        << _name << ":" << pyre::journal::newline
        << "      iteration: " << _iteration << pyre::journal::newline
        << "        workers: " << _workers << pyre::journal::newline
        << "  chains/worker: " << annealer.samples() << pyre::journal::newline
        << "   total chains: " << _workers * annealer.samples() << pyre::journal::newline;

    // all done
    return *this;
}

//// @par Main functionality
/// Read beta and scaling factor from BetaStatistics file
/// @param [in] file File name of the sample statistics file from a previous run
/// @param [in] iteration The iteration number to be restarted
/// @param [in, out] beta The annealing temperature to be read in
/// @param [in, out] scaling The random walk scaling factor to be read in
/// @note
/// This function is privately used by the "restart" function of the same class
altar::bayesian::AnnealingMethod &
altar::bayesian::AnnealingMethod::
_readBetastatistics(std::string file, const size_t iteration, TYPE & beta, TYPE & scaling)
{
    std::ifstream fin(file.c_str());
    if (!fin) std::cerr<<"Cannot open "<<file<<std::endl;
    std::string line;
    std::stringstream ss; 
    std::string word;
    size_t currentIteration;
    std::string dummy;
    while (getline(fin, line))
    {   
        ss.str(line);
        if (ss>>currentIteration, currentIteration==iteration)
        {   
            ss>>scaling>>dummy>>beta;
            return *this;
        }   
        ss.clear();
    }   
    fin.close();
    return *this;
}


// end of file
