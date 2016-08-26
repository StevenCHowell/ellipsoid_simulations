#include <iostream>
#include <ellipsoid.h>


////////////////////////////////////////////////////////////////////////
// scale and truncate energy values
////////////////////////////////////////////////////////////////////////
void process_energy(Eigen::VectorXf & energy, const float Umin)
{
    bool flag = false;
    energy /= Umin;
    for (int i=energy.size()-2; i>=0; --i)
    {
        if (energy(i)>=energy(i+1)) flag = true;
        else if (flag) energy(i) = energy(i+1);
    }
}


////////////////////////////////////////////////////////////////////////
// main test
////////////////////////////////////////////////////////////////////////
void reproduce_fig_8()
{
    // force field parameters
    const float ep = 25.543;
    const float sig = 2.921;

    // the repulsive shape matrix elements for A1 and B1
    const float a11=7.024;
    const float a12=4.006;
    const float a13=1.582;
    const float b11=4.846;
    const float b12=2.841;
    const float b13=1.511;
    // the attractive shape matrix elements for A2 and B2
    const float a21=6.528;
    const float a22=3.829;
    const float a23=1.702;
    const float b21=4.267;
    const float b22=2.405;
    const float b23=1.394;

	Eigen::Vector3f r; r<<0,0,0;//particle 1 is at the origin
	Eigen::VectorXf sx_coord = Eigen::ArrayXf::LinSpaced(300,3.0,16.0); //x-coordinate of particle 2

    /*
    // 1-1
	Eigen::Vector3f axis_lengths_r1; axis_lengths_r1<<a11,a12,a13;
	Eigen::Vector3f axis_lengths_r2; axis_lengths_r2<<a21,a22,a23;
	Eigen::Vector3f axis_lengths_s1; axis_lengths_s1<<b11,b12,b13;
	Eigen::Vector3f axis_lengths_s2; axis_lengths_s2<<b21,b22,b23;
    // 1-2
	Eigen::Vector3f axis_lengths_r1; axis_lengths_r1<<a11,a12,a13;
	Eigen::Vector3f axis_lengths_r2; axis_lengths_r2<<a21,a22,a23;
	Eigen::Vector3f axis_lengths_s1; axis_lengths_s1<<b12,b11,b13;
	Eigen::Vector3f axis_lengths_s2; axis_lengths_s2<<b22,b21,b23;
    // 1-3
	Eigen::Vector3f axis_lengths_r1; axis_lengths_r1<<a11,a12,a13;
	Eigen::Vector3f axis_lengths_r2; axis_lengths_r2<<a21,a22,a23;
	Eigen::Vector3f axis_lengths_s1; axis_lengths_s1<<b13,b11,b12;
	Eigen::Vector3f axis_lengths_s2; axis_lengths_s2<<b23,b21,b22;
    // 2-1
	Eigen::Vector3f axis_lengths_r1; axis_lengths_r1<<a12,a11,a13;
	Eigen::Vector3f axis_lengths_r2; axis_lengths_r2<<a22,a21,a23;
	Eigen::Vector3f axis_lengths_s1; axis_lengths_s1<<b11,b12,b13;
	Eigen::Vector3f axis_lengths_s2; axis_lengths_s2<<b21,b22,b23;
    // 2-2
	Eigen::Vector3f axis_lengths_r1; axis_lengths_r1<<a12,a11,a13;
	Eigen::Vector3f axis_lengths_r2; axis_lengths_r2<<a22,a21,a23;
	Eigen::Vector3f axis_lengths_s1; axis_lengths_s1<<b12,b11,b13;
	Eigen::Vector3f axis_lengths_s2; axis_lengths_s2<<b22,b21,b23;
    */
    // 2-3
	Eigen::Vector3f axis_lengths_r1; axis_lengths_r1<<a12,a11,a13;
	Eigen::Vector3f axis_lengths_r2; axis_lengths_r2<<a22,a21,a23;
	Eigen::Vector3f axis_lengths_s1; axis_lengths_s1<<b13,b11,b12;
	Eigen::Vector3f axis_lengths_s2; axis_lengths_s2<<b23,b21,b22;
    /*
    // 3-1
	Eigen::Vector3f axis_lengths_r1; axis_lengths_r1<<a13,a11,a12;
	Eigen::Vector3f axis_lengths_r2; axis_lengths_r2<<a23,a21,a22;
	Eigen::Vector3f axis_lengths_s1; axis_lengths_s1<<b11,b12,b13;
	Eigen::Vector3f axis_lengths_s2; axis_lengths_s2<<b21,b22,b23;
    // 3-2
	Eigen::Vector3f axis_lengths_r1; axis_lengths_r1<<a13,a11,a12;
	Eigen::Vector3f axis_lengths_r2; axis_lengths_r2<<a23,a21,a22;
	Eigen::Vector3f axis_lengths_s1; axis_lengths_s1<<b12,b11,b13;
	Eigen::Vector3f axis_lengths_s2; axis_lengths_s2<<b22,b21,b23;
    // 3-3
	Eigen::Vector3f axis_lengths_r1; axis_lengths_r1<<a13,a11,a12;
	Eigen::Vector3f axis_lengths_r2; axis_lengths_r2<<a23,a21,a22;
	Eigen::Vector3f axis_lengths_s1; axis_lengths_s1<<b13,b11,b12;
	Eigen::Vector3f axis_lengths_s2; axis_lengths_s2<<b23,b21,b22;
    */

	const int len = sx_coord.size();

	Eigen::VectorXf energy(len);
	Eigen::VectorXf virial(len);
	Eigen::Vector3f force, torque;

	const float boxsize = 1000;

    ellipsoid::Ellipsoid e1(0,0,0,a12,a11,a13,a22,a21,a23,ep,sig);
    ellipsoid::Ellipsoid e2(0,0,0,b13,b11,b12,b23,b21,b22,ep,sig);
    //ellipsoid::Ellipsoid e1(r,axis_lengths_r1, axis_lengths_r2, ep, sig);
    //ellipsoid::Ellipsoid e2(r,axis_lengths_s1, axis_lengths_s2, ep, sig);

    Eigen::AngleAxisf rot(0.0, Eigen::Vector3f(0.0,1.0,0.0));
	for (int i=0; i<len; ++i)
	{
		Eigen::Vector3f s; s<<sx_coord(i),0,0;
        e1.set_Avec(r);
        e2.set_Avec(s);
        e1.rotate_Amat(rot.matrix());
        e2.rotate_Amat(rot.matrix());
		e1.energy_pbc(e2, boxsize, energy(i));
	}

    process_energy(energy, 25.0);
    for (int i=0; i<len; ++i) std::cout<<sx_coord(i)<<" "<<energy(i)<<std::endl;
};

////////////////////////////////////////////////////////////////////////
// main test
////////////////////////////////////////////////////////////////////////
int main()
{
    reproduce_fig_8();
}
