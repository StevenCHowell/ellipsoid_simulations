#include "ellipsoid.h"

#include <math.h>
#include <iostream>


/// @par Main functionality
/// Calculate the FAB value<br>
/// @param Amat is the matrix for ellipsoid 1 in Eq. 1 
/// @param Bmat is the matrix for ellipsoid 1 in Eq. 1
/// @note all of the above quantities are defined as matrices, not arrays
float
FAB::
calc_FAB(const Eigen::Vector3f & Avec, const Eigen::Vector3f & Bvec, const Eigen::Matrix3f & Amat, const Eigen::Matrix3f & Bmat, Eigen::Vector3f & x_of_lam_new, float & lam_new)
{
    FAB::FAB_params p;
    p.p_Avec = &Avec;
    p.p_Bvec = &Bvec;
    p.p_Amat = &Amat;
    p.p_Bmat = &Bmat;

    int status;
    int iter = 0, max_iter = 100;
    const gsl_root_fsolver_type *T;
    gsl_root_fsolver *s;
    double lam = 0.5;
    double lam_low = 0.0, lam_high = 1.0;
    gsl_function F;

    F.function = &(FAB::func_FAB);
    F.params = &p;

    T = gsl_root_fsolver_brent;
    s = gsl_root_fsolver_alloc (T);
    gsl_root_fsolver_set (s, &F, lam_low, lam_high);

    //printf ("using %s method\n", gsl_root_fsolver_name (s));
    //printf ("%5s [%9s, %9s] %9s %10s %9s\n", "iter", "lower", "upper", "root", "err", "err(est)");

    do
    {
        iter++;
        status = gsl_root_fsolver_iterate (s);
        lam = gsl_root_fsolver_root (s);
        lam_low = gsl_root_fsolver_x_lower (s);
        lam_high = gsl_root_fsolver_x_upper (s);
        status = gsl_root_test_interval (lam_low, lam_high, 0, 0.001);

        //if (status == GSL_SUCCESS) printf ("Converged:\n");
        //printf ("%5d [%.7f, %.7f] %.7f %.7f\n", iter, lam_low, lam_high, lam, lam_high - lam_low);
    }
    while (status == GSL_CONTINUE && iter < max_iter);

    gsl_root_fsolver_free (s);

    x_of_lam_new = p.x_of_lam;
    lam_new = lam;

    return p.Ascr;
}

/// @par Main functionality
/// function to calculate FAB (used for gsl minimization) <br>
double 
FAB::
func_FAB(double lam, void * params)
{
    FAB_params & p = *static_cast<FAB_params*>(params);
    const Eigen::Vector3f & Avec = *(p.p_Avec);
    const Eigen::Matrix3f & Amat = *(p.p_Amat);
    const Eigen::Vector3f & Bvec = *(p.p_Bvec);
    const Eigen::Matrix3f & Bmat = *(p.p_Bmat);
    p.x_of_lam = ((lam*Amat+(1.0-lam)*Bmat).inverse()) * (lam*Amat*Avec + (1.0-lam)*Bmat*Bvec);
    const float Ascr = ((p.x_of_lam-Avec).transpose())*Amat*(p.x_of_lam-Avec);
    const float Bscr = ((p.x_of_lam-Bvec).transpose())*Bmat*(p.x_of_lam-Bvec);
    p.Ascr = Ascr;
    return Ascr-Bscr;
}



/// @par Main functionality
/// calculate energy based on periodic boundary conditions <br>
/// calculates the interaction energy due to two ellipsoids using the method from Paramonov and Yaliraki, 2005<br>
/// @param box_length is the simulation box length
/// @return the energy
bool
ellipsoid::Ellipsoid::
energy_pbc(const ellipsoid::Ellipsoid & e, const float box_length, float & energy) const
{
    Eigen::Vector3f Bvec_pbc = e.Avec()-_Avec;
    for (int i=0; i<3; ++i) Bvec_pbc[i] -= box_length*(round(Bvec_pbc[i]/box_length)) ;
    Bvec_pbc += _Avec;

	const float distR = (Bvec_pbc-_Avec).norm(); // the distance between the ellipsoids given the periodic boundaries, 'R' in the paper
	
	//std::cout<<"Avec: "<<_Avec<<std::endl<<"Bvec_pbc: "<<Bvec_pbc<<std::endl<<" distR= "<<distR<<std::endl;
	// std::cout<<"distR="<<distR<<std::endl;

    Eigen::Vector3f xc1;
	float lam1;
	const float FAB1 = FAB::calc_FAB(_Avec,Bvec_pbc,_Amat1,e.Amat1(),xc1,lam1);
    const float denominator1 = distR-distR*pow(FAB1,-0.5)+_sig;
    if (denominator1<=_machine_precision_float)
    {
        energy = _machine_inf_float;
        return false;
    }

    Eigen::Vector3f xc2;
	float lam2;
	const float FAB2 = FAB::calc_FAB(_Avec,Bvec_pbc,_Amat2,e.Amat2(),xc2,lam2);
    const float denominator2 = distR-distR*pow(FAB2,-0.5)+_sig;
    if (denominator2<=_machine_precision_float)
    {
        energy = _machine_inf_float;
        return false;
    }

	energy = 4*_ep*(pow(_sig/denominator1,12.0)-pow(_sig/denominator2,6.0)); // Eq. 44
	//energy = 4*_ep*(pow(_sig/(distR-distR*pow(FAB1,-0.5)+_sig),12.0)-pow(_sig/(distR-distR*pow(FAB2,-0.5)+_sig),6.0)); // Eq. 44

	return true;
}

/// @par Main functionality
/// destructor <br>
ellipsoid::Ellipsoid::
~Ellipsoid()
{
}
