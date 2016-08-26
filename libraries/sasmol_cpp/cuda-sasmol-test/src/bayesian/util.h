#if !defined(altar_bayesian_util_h)
#define altar_bayesian_util_h

// macros
#include <altar/utils/common.h>

#include <altar/utils/util.h>

#include <stdio.h>
#include <string.h>
#include <string>
#include <sys/time.h>
#include <sstream>
#include <sys/stat.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_eigen.h>
#ifdef USE_DOUBLE 
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#else
#include <gsl/gsl_vector_float.h>
#include <gsl/gsl_matrix_float.h>
#endif


// place everything in the local namespace
namespace altar {
    namespace bayesian {
		namespace util {

			// condition a matrix to be positive definate
			void conditionMatrix(gsl_matrix_TYPE * M);

		} //namespace util
	} // namespace bayesian
} // namespace altar

#endif
