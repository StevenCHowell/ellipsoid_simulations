#include "util.h"

////////////////////////////////////////////////////////////
// condition a matrix to be positive definate
////////////////////////////////////////////////////////////
void
altar::bayesian::util::
conditionMatrix(gsl_matrix_TYPE * M)
{
	size_t m = M->size1;
	TYPE eval_ratio_min = 0.001;
     
    gsl_vector_TYPE *eval = gsl_vector_TYPE_alloc (m);
    gsl_matrix_TYPE *evec = gsl_matrix_TYPE_alloc (m, m);     
    gsl_eigen_symmv_workspace * w = gsl_eigen_symmv_alloc (m);
#ifdef USE_DOUBLE
    gsl_eigen_symmv (M, eval, evec, w);
    gsl_eigen_symmv_free (w);
    gsl_eigen_symmv_sort (eval, evec,  GSL_EIGEN_SORT_ABS_ASC);
#else
	gsl_matrix * M_double = altar::utils::gsl_matrix_double_from_float(M);
	gsl_vector * eval_double = altar::utils::gsl_vector_double_from_float(eval);
	gsl_matrix * evec_double = altar::utils::gsl_matrix_double_from_float(evec);
    gsl_eigen_symmv (M_double, eval_double, evec_double, w);
    gsl_eigen_symmv_free (w);
    gsl_eigen_symmv_sort (eval_double, evec_double,  GSL_EIGEN_SORT_ABS_ASC);
	altar::utils::gsl_matrix_set_float_from_double(M, M_double);
	altar::utils::gsl_vector_set_float_from_double(eval, eval_double);
	altar::utils::gsl_matrix_set_float_from_double(evec, evec_double);
	gsl_matrix_free(M_double);
	gsl_vector_free(eval_double);
	gsl_matrix_free(evec_double);
#endif

	gsl_matrix_TYPE *evecT = gsl_matrix_TYPE_alloc(m,m);
	gsl_matrix_TYPE_transpose_memcpy(evecT, evec);
	gsl_matrix_TYPE *diagM = gsl_matrix_TYPE_calloc(m,m);
	TYPE eval_min = eval_ratio_min*gsl_vector_TYPE_get(eval,m-1);
    TYPE eval_i;
//printf("eigenvalues:\n");
    for (size_t i = 0; i < m; i++)
    {
    	eval_i  = gsl_vector_TYPE_get (eval, i);
//printf ("%f ", eval_i);
		if (eval_i<eval_min) gsl_matrix_TYPE_set(diagM, i, i, eval_min);
		else gsl_matrix_TYPE_set(diagM, i, i, eval_i);
    }
//printf("\n");
//mprint(diagM, "diagM:", 15, 15);
	gsl_matrix_TYPE *tmp = gsl_matrix_TYPE_alloc(m,m);
	gsl_blas_TYPEgemm(CblasNoTrans, CblasNoTrans, 1.0, diagM, evecT, 0.0, tmp);
	gsl_blas_TYPEgemm(CblasNoTrans, CblasNoTrans, 1.0, evec, tmp, 0.0, M);

	gsl_matrix_TYPE_transpose_memcpy(tmp, M);
	gsl_matrix_TYPE_add(M,tmp);
	gsl_matrix_TYPE_scale(M, 0.5);
   
    gsl_vector_TYPE_free (eval);
    gsl_matrix_TYPE_free (evec);
    gsl_matrix_TYPE_free (evecT);
    gsl_matrix_TYPE_free (diagM);
    gsl_matrix_TYPE_free (tmp);
}


