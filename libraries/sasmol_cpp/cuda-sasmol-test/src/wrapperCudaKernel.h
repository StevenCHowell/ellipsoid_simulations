//
// Hailiang Zhang
// NIST & UTK
//

#ifndef H_WAPPERCUDAKERNEL
#define H_WAPPERCUDAKERNEL

#include "cudaUtil.h"
#include <curand_kernel.h>

namespace sascalc
{
    namespace cuda
    {
        void wrapper_cudaMC_setupCurandState(curandState * const state, const int seed, const int Npoints);

        void wrapper_cudaMC_fill(const float * const gpu_coor, const float * const gpu_radii, float * const gpu_points,
                                const int Natoms, const int Npoints,
                                const float x0, const float y0, const float z0, const float xl, const float yl, const float zl,
                                curandState * const state);

        void wrapper_cudaGV_calcIq(const float * const gpu_coor, const float * const gpu_b, const float * const gpu_gv, const float * const gpu_q, float * const gpu_Is, const int Natoms, const int Nq, const int NGV);

        void wrapper_cudaDebye_calcIq(const float * const gpu_coor, const float * const gpu_b, const float * const gpu_q, float * const gpu_Ia, const int Natoms, const int Nq);
    }
}

#endif
