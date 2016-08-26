#ifndef H_ELLIPSOID
#define H_ELLIPSOID

#include <vector>
#include <cmath>
#include <Eigen/Dense>

typedef Eigen::Array<float,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor> EigenArrayXXfRow;
typedef Eigen::Array<float,1,3> EigenArray13fRow;

#define A1 1.0
#define A1 1.0
#define N1d 5 
#define RADIUS 1.0
#define RSQCUT (RADIUS*RADIUS*4)
#define RBEADTOBEAD (RADIUS*3.1934)
#define RCENTERTOBEAD (RADIUS*1.8437)

#define RHH RSQCUT // cut off radius
#define RHL1 RSQCUT // cut off radius
#define RHL2 RSQCUT // cut off radius
#define RL1L1 RSQCUT // cut off radius
#define RL1L2 RSQCUT // cut off radius
#define RL2L2 RSQCUT // cut off radius

// globals aren't good but I will move it later
const EigenArray13fRow VECTOR_CENTER_TO_H(RCENTERTOBEAD, 0.0, 0.0);
const EigenArray13fRow VECTOR_CENTER_TO_L1(-RCENTERTOBEAD/2.0, RCENTERTOBEAD*1.732/2.0, 0.0);
const EigenArray13fRow VECTOR_CENTER_TO_L2(-RCENTERTOBEAD/2.0, -RCENTERTOBEAD*1.732/2.0, 0.0);
const EigenArray13fRow U0_H(1.,0.,0.);
const EigenArray13fRow U1_H(0.,1.,0.);
const EigenArray13fRow U2_H(0.,0.,1.);
const EigenArray13fRow U0_L1(-0.5,1.7320508075688772,0.);
const EigenArray13fRow U1_L1(-1.7320508075688772,-0.5,0.);
const EigenArray13fRow U2_L1(0.,0.,1.);
const EigenArray13fRow U0_L2(-0.5,-1.7320508075688772,0.);
const EigenArray13fRow U1_L2(1.7320508075688772,-0.5,0.);
const EigenArray13fRow U2_L2(0.,0.,1.);

// ellipsoid paramters
#define RCENTERTOH RCENTERTOBEAD
#define RCENTERTOL1 RCENTERTOBEAD
#define RCENTERTOL2 RCENTERTOBEAD
const Eigen::Vector3f axis_lengths_H(1,1,1);
const Eigen::Vector3f axis_lengths_L1(1,1,1);
const Eigen::Vector3f axis_lengths_L2(1,1,1);

const float sigma = 1.0;
const float epsilon = 1000.0;

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
               const float sig, const float ep, const float boxsize);

#endif
