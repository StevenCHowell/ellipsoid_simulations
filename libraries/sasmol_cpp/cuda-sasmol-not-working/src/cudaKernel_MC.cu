//
// Hailiang Zhang
// NIST & UTK
//

#include "cudaKernel_MC.h"


/// @par Main functionality
/// Check whether a sample is a bad move for the prior
/// @param [in] gM_tmp the candiate @b M matrix (Dimension: parameters*samples, samples is the leading index)
__global__ void setupCurandKernel(curandState * const state, const int seed, const int n)
{
    int id = threadIdx.x + blockIdx.x*blockDim.x;
    if(id<n) curand_init(seed, id, 0, &state[id]);
}

/// @par Main functionality
/// Check whether a sample is a bad move for the prior
/// @param [in] gM_tmp the candiate @b M matrix (Dimension: parameters*samples, samples is the leading index)
__global__ void cudaMC(const float * const gpu_coor, const float * const gpu_radii, float * const gpu_points,
                       const int Natoms, const int Npoints,
                       const float x0, const float y0, const float z0, const float xl, const float yl, const float zl,
                       curandState * const state)
{
	// get the thread id
    const int idx = blockIdx.x*blockDim.x + threadIdx.x;
    if (idx >= Npoints) return;

    //if (idx<10) printf("idx: %d,  Natoms: %d Npoints: %d,  x0: %f y0: %f z0: %f,  xl: %f yl: %f zl: %f,  \n",idx,Natoms,Npoints,x0,y0,z0,xl,yl,zl);

    curandState localState = state[idx];

    float x,y,z;
    float r;
    bool flag_found = false;
    while (!flag_found)
    {
        // get a new point
        x = x0 + curand_uniform(&localState)*xl;
        y = y0 + curand_uniform(&localState)*yl;
        z = z0 + curand_uniform(&localState)*zl;
        //if (idx<10) printf("idx: %d,  Natoms: %d Npoints: %d,  x0: %f y0: %f z0: %f,  xl: %f yl: %f zl: %f,  x: %f y: %f z:%f\n",idx,Natoms,Npoints,x0,y0,z0,xl,yl,zl,x,y,z);
        // check whether this points is inside the molecule or not
        for (int i=0; i<Natoms; ++i)
        {
            r = sqrt(powf((x-gpu_coor[i]),2.0)+powf((y-gpu_coor[Natoms+i]),2.0)+powf((z-gpu_coor[2*Natoms+i]),2.0));
            if (r<gpu_radii[i])
            {
                flag_found = true;
                gpu_points[idx] = x;
                gpu_points[Npoints+idx] = y;
                gpu_points[2*Npoints+idx] = z;
                break;
            }
        }
    }
}
