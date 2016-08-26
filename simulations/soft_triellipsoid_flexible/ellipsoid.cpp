#include "ellipsoid.h"

#include <math.h>
#include <iostream>
#include <Eigen/Dense>


// my sign template
template <class T> inline int sign(T val) {return (T(0)<val)-(T(0)>val);}


////////////////////////////////////////////////////////////////////////////////////////////////
// calculates the separation vector that points from r1 to r2
// using the minimum image convention of an orthorhombic box
// boxsize can be a single scalar or a 3D vector
////////////////////////////////////////////////////////////////////////////////////////////////
Eigen::Vector3f dx_pbc(const Eigen::Vector3f & r1, const Eigen::Vector3f & r2, const float boxsize)
{
    Eigen::Vector3f dx;
    for (int i=0; i<3; ++i)
    {
        dx[i] = r1[i] - floor(r1[i]/boxsize)*boxsize;
        dx[i] = r2[i] - dx[i];
        dx[i] -= round(dx[i]/boxsize)*boxsize;
    }
    return dx;
}


////////////////////////////////////////////////////////////////////////////////////////////////
// returns the square of the distance between 2 locations
// whose separation is calculated in dx_pbc
////////////////////////////////////////////////////////////////////////////////////////////////
float dist2_pbc(const Eigen::Vector3f & r1, const Eigen::Vector3f & r2, const float boxsize)
{
    Eigen::Vector3f dx = dx_pbc(r1, r2, boxsize);
    return dx[0]*dx[0] + dx[1]*dx[1] + dx[2]*dx[2];
}

////////////////////////////////////////////////////////////////////////////////////////////////
// calculates the matrix which describes the shape of an ellipsoid according to Eq. 1
// axis_lengths is a 3-element array containing the length of the 3 semiaxes
// u0,1,2 are 3 orthonormal unit vectors which describe the orientation of the 3 axes. 
////////////////////////////////////////////////////////////////////////////////////////////////
Eigen::Matrix3f ellipsoid_matrix(const Eigen::Vector3f & axis_lengths, const Eigen::Vector3f & u0, const Eigen::Vector3f & u1, const Eigen::Vector3f & u2)
{
	const float a0 = float(axis_lengths[0]);
	const float a1 = float(axis_lengths[1]);
	const float a2 = float(axis_lengths[2]);

	Eigen::Matrix3f output = u0*u0.transpose()/(a0*a0) + u1*u1.transpose()/(a1*a1) + u2*u2.transpose()/(a2*a2);

	return output;
}

////////////////////////////////////////////////////////////////////////////////////////////////
// calculates x(lambda) according to Eq. 5 in the form of a column vector
////////////////////////////////////////////////////////////////////////////////////////////////
Eigen::Vector3f calc_x_of_lam(const Eigen::Vector3f & rvect, const Eigen::Matrix3f & Amat, const Eigen::Vector3f & svect, const Eigen::Matrix3f & Bmat, const float lam)
{

	Eigen::Vector3f output = ((lam*Amat+(1.0-lam)*Bmat).inverse()) * (lam*Amat*rvect + (1.0-lam)*Bmat*svect);
	// std::cout<<"x_of_lam="<<output<<std::endl;

	return output;
}


////////////////////////////////////////////////////////////////////////////////////////////////
// rvect is the com position of ellipsoid 1
// svect is the com position of ellipsoid 2
// Amat is the matrix for ellipsoid 1 in Eq. 1 
// Bmat is the matrix for ellipsoid 1 in Eq. 1
//  all of the above quantities are defined as matrices, not arrays
////////////////////////////////////////////////////////////////////////////////////////////////
float calc_F_AB(const Eigen::Vector3f & rvect, const Eigen::Vector3f & svect, const Eigen::Matrix3f & Amat, const Eigen::Matrix3f & Bmat, Eigen::Vector3f & x_of_lam_new)
{

	const float tol=0.001; // tolerance for finding lambda, which ranges between 0 and 1
	const int Nmax=50; // maximum number of steps for root-finding
	
	Eigen::Vector3f x_of_lam0=calc_x_of_lam(rvect,Amat,svect,Bmat,0.0); // lambda=0
	Eigen::Vector3f x_of_lam1=calc_x_of_lam(rvect,Amat,svect,Bmat,1.0); // lambda=1
	
	const float fa = float(((x_of_lam0-rvect).transpose())*Amat*(x_of_lam0-rvect)) - float(((x_of_lam0-svect).transpose())*Bmat*(x_of_lam0-svect)); // f(a), when lambda=0 and f=A-B in Eq. 9. with scripted A and B before Eq. 1
	const float fb = float(((x_of_lam1-rvect).transpose())*Amat*(x_of_lam1-rvect)) - float(((x_of_lam1-svect).transpose())*Bmat*(x_of_lam1-svect)); // f(b)

	int Nsteps = 1;
	bool converged = false;
	
	float a=0.0, b=1.0; // the lower and upper bounds for lambda
	
    float Ascr, Bscr, fc, c;

	while(converged==false) // use bisection method for finding 
    {
	
		// std::cout<<"N="<<Nsteps<<std::endl;
	
		if(Nsteps==Nmax)
        {
            std::cout<<"maximum number of steps reached in root finding"<<std::endl;
            exit(0);
        }
		
		c = 0.5*(b+a);
		
	    // std::cout<<"lam="<<c<<std::endl;
		
		x_of_lam_new = calc_x_of_lam(rvect,Amat,svect,Bmat,c);
		
		Ascr = float(((x_of_lam_new-rvect).transpose())*Amat*(x_of_lam_new-rvect));
		Bscr = float(((x_of_lam_new-svect).transpose())*Bmat*(x_of_lam_new-svect));
		// the scripted A and B defined before Eq. 1
		
		fc = Ascr - Bscr; // f(c)
		
		if (fc==0 || 0.5*(b-a)<tol) converged = true;
		else
        {
			if (sign(fc) == sign(fa)) a = c;
			else b = c;
        }

		Nsteps++;
    }

	//output=Ascr; // from Eq. 11

	return Ascr;
}

////////////////////////////////////////////////////////////////////////////////////////////////
// calculates the interaction energy due to two ellipsoids using the method from Paramonov and Yaliraki, 2005

// r is the com position of ellipsoid 1, which is a numpy array
// s is the com position of ellipsoid 2, which is a numpy array

// axis_lengths_r is a 3-element array containing the lengths of the 3 semiaxes of ellipsoid 1
// axis_lengths_s is a 3-element array containing the lengths of the 3 semiaxes of ellipsoid 2

// u0,1,2 are 3 orthonormal unit vectors which describe the 3 axes of ellipsoid 1
// v0,1,2 are 3 orthonormal unit vectors which describe the 3 axes of ellipsoid 2
// sig and ep are the Lennard-Jones parameters
// boxsize is the simulation box size
////////////////////////////////////////////////////////////////////////////////////////////////
float pot_energy(const Eigen::Vector3f & r, const Eigen::Vector3f & axis_lengths_r, const Eigen::Vector3f & u0, const Eigen::Vector3f & u1, const Eigen::Vector3f & u2,
               const Eigen::Vector3f & s, const Eigen::Vector3f & axis_lengths_s, const Eigen::Vector3f & v0, const Eigen::Vector3f & v1, const Eigen::Vector3f & v2,
               const float sig, const float ep, const float boxsize)
{

	const float distR = sqrt(dist2_pbc(r,s,boxsize)); // the distance between the ellipsoids given the periodic boundaries, 'R' in the paper
	
	// std::cout<<"distR="<<distR<<std::endl;

	const Eigen::Matrix3f Amat = ellipsoid_matrix(axis_lengths_r,u0,u1,u2);
	const Eigen::Matrix3f Bmat = ellipsoid_matrix(axis_lengths_s,v0,v1,v2);
	
	// convert arrays to column vectors, which are numpy matrices
	
	const Eigen::Vector3f rvect = r;
	
	const Eigen::Vector3f svect = s;
	
    Eigen::Vector3f xc;
	const float F_AB = calc_F_AB(rvect,svect,Amat,Bmat,xc);

	const Eigen::Vector3f xa = rvect+(xc-rvect)/sqrt(F_AB);
	const Eigen::Vector3f xb = svect+(xc-svect)/sqrt(F_AB);

	// std::cout<<"xa="<<xa<<std::endl;
	// std::cout<<"xb="<<xb<<std::endl;
	
	float term=(sig/(distR-distR/sqrt(F_AB)+sig));
	term=term*term*term;
	term=term*term; // becomes to the 6th power
	
	const float U=4*ep*(term*term-term); // Eq. 43

	return U;
}

/*
////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////
// main entrance
////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////
int main()
{
    Eigen::Vector3f u0(1,0,0);
    Eigen::Vector3f u1(0,1,0);
    Eigen::Vector3f u2(0,0,1);

    Eigen::Vector3f v0(1,0,0);
    Eigen::Vector3f v1(0,1,0);
    Eigen::Vector3f v2(0,0,1);

    Eigen::Vector3f axis_lengths_r(1,1,1);
    Eigen::Vector3f axis_lengths_s(1,1,1);

    Eigen::VectorXf sx_coord(170);
    float f=1.5;
    for (int i=0; f<10.; ++i)
    {
        sx_coord(i)=f;
        f+=0.05;
    }
    Eigen::VectorXf energy;
    energy.setZero(sx_coord.size());

for (int count=0; count<100; ++count)
    for (int i=0; i<sx_coord.size(); ++i)
    {
        Eigen::Vector3f r(0,0,0);
        Eigen::Vector3f s(sx_coord[i],0,0);

        float sig=1.0;
        float ep=1.0;
        float boxsize=1000;

        energy[i] = pot_energy(r,axis_lengths_r,u0,u1,u2,s,axis_lengths_s,v0,v1,v2,sig,ep,boxsize);
    }
    for (int i=0; i<sx_coord.size(); ++i)
        std::cout<<sx_coord[i]<<" "<<energy[i]<<std::endl;

    return 0;
}
*/
