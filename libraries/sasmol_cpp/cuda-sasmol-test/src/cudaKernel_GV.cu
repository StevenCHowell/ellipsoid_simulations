//
// Hailiang Zhang
// NIST & UTK
//

#include "cudaKernel_GV.h"


/// @par Main functionality
/// Calculate I(q) based on GV method
/// @param [in] gpu_coor the XYX coordinates where the leading dimension is Natoms and the slow dimension is X->Y->Z
/// @param [in] gpu_gv the gv data where the leading dimension is NGV and the slow dimension is X->Y->Z
/// @param [in] gpu_Ireal the Ireal data where the leading dimension is Nq and the slow dimension is NGV
__global__ void
sascalc::cuda::
cudaKernel_GV_calc(const float * const gpu_coor, const float * const gpu_b, const float * const gpu_gv, const float * const gpu_q, float * const gpu_I, const int Natoms, const int Nq, const int NGV)
{
	// get the thread id
    //int idx = blockIdx.x*blockDim.x + threadIdx.x;

    //printf("blockIdx.x: %d blockDim.x: %d threadIdx.x: %d idx: %d Natoms %d Nq %d NGV %d\n",blockIdx.x,blockDim.x,threadIdx.x,idx,Natoms,Nq,NGV);

    // get q values
    const int iq = threadIdx.x;
    const int igv = blockIdx.x;
    const float qmag = gpu_q[iq];
    const float qx = qmag * gpu_gv[igv];
    const float qy = qmag * gpu_gv[NGV + igv];
    const float qz = qmag * gpu_gv[NGV*2 + igv];

    // locals
    int iatom;
    float x, y, z, b;
    float q_dot_r;
    float sine, cose;
    
    // summation
    float Ireal = 0.0;
    float Iimag = 0.0;
    for (iatom=0; iatom<Natoms; ++iatom)
    {
        x = gpu_coor[iatom];
        y = gpu_coor[Natoms + iatom];
        z = gpu_coor[Natoms*2 + iatom];
        b = gpu_b[iatom];
        q_dot_r = qx*x+qy*y+qz*z;
        sincosf(q_dot_r, &sine, &cose);
        Ireal += b*cose;
        Iimag += b*sine;
    }

    // save
    gpu_I[igv*Nq + iq] = Ireal*Ireal + Iimag*Iimag;
}
