//
// Hailiang Zhang
// NIST & UTK
//

#include "cudaKernel_Debye.h"


/// @par Main functionality
/// Calculate I(q) based on Debye summation
/// @param [in] gpu_coor the XYX coordinates where the leading dimension is Natoms and the slow dimension is X->Y->Z
/// @param [in] gpu_Ia the partial intensity where the leading dimension is Natoms and the slow dimension is Nq
__global__ void
sascalc::cuda::
cudaKernel_Debye_calc(const float * const gpu_coor, const float * const gpu_b, const float * const gpu_q, float * const gpu_Ia, const int Natoms, const int Nq)
{
	// get the thread id and return if necessary
    const int idx = blockIdx.x*blockDim.x + threadIdx.x;
    if (idx>=Natoms) return;

    // shared memory
    __shared__ float shared_x[BLOCKDIM];
    __shared__ float shared_y[BLOCKDIM];
    __shared__ float shared_z[BLOCKDIM];
    __shared__ float shared_b[BLOCKDIM];

    // coordinate/b for this thread
    const float xi = gpu_coor[idx];
    const float yi = gpu_coor[Natoms+idx];
    const float zi = gpu_coor[Natoms*2+idx];
    const float bi = gpu_b[idx];

    // local variables
    int idx_atom,idx_block, idx_q;
    float r,q,qr,bj;

                for (idx_q=0; idx_q<Nq; ++idx_q) gpu_Ia[idx_q*Natoms + idx] = 0.0;

    // loop over number of blocks
    int idx0_atom = threadIdx.x; // the first atom_j to start with
    //for (idx_block=blockIdx.x; idx_block<((Natoms-1)/blockDim.x)+1; ++idx_block)
    for (idx_block=blockIdx.x; idx_block<gridDim.x; ++idx_block)
    {
        // load coordinates and bs into shared memory
        __syncthreads();
        idx_atom = idx_block*blockDim.x + threadIdx.x;
        shared_x[threadIdx.x] = gpu_coor[idx_atom];
        shared_y[threadIdx.x] = gpu_coor[Natoms + idx_atom];
        shared_z[threadIdx.x] = gpu_coor[Natoms*2 + idx_atom];
        shared_b[threadIdx.x] = gpu_b[idx_atom];
        __syncthreads();
        for (idx_atom=idx0_atom; idx_atom<blockDim.x; ++idx_atom)
        {
            if (idx_block*blockDim.x + idx_atom >= Natoms) break;
            r = sqrt ( powf(xi-shared_x[idx_atom], 2.0) + powf(yi-shared_y[idx_atom], 2.0) + powf(zi-shared_z[idx_atom], 2.0) );
            bj = shared_b[idx_atom];
            for (idx_q=0; idx_q<Nq; ++idx_q)
            {
                q = gpu_q[idx_q];
                qr = q*r;
                if (qr==0.0) gpu_Ia[idx_q*Natoms + idx] += bi*bj;
                else gpu_Ia[idx_q*Natoms + idx] += bi*bj*sinf(q*r)/(q*r);
            }
         }
         idx0_atom = 0;
    }
}
