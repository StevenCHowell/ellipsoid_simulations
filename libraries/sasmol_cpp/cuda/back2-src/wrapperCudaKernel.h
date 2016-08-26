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
        void wrapper_cudaMC_setupGrid(const float * const gpu_coor, const float * const gpu_radii, int * const gpu_grid, const int Natoms,
                                      const float x0, const float y0, const float z0, const float delta, const int Nx, const int Ny, const int Nz);
        void wrapper_cudaMC_fill(const float * const gpu_coor, const float * const gpu_radii, const int * const gpu_grid, float * const gpu_points,
                                const int Natoms, const int Npoints, int * const gpu_Npoints_tried,
                                const float x0, const float y0, const float z0, const float delta, const int Nx, const int Ny, const int Nz,
                                curandState * const state);

        void wrapper_cudaGV_calcIq(const float * const gpu_coor, const float * const gpu_b, const float * const gpu_gv, const float * const gpu_q, float * const gpu_Is, const int Natoms, const int Nq, const int Ngv);
        void wrapper_cudaGV_calcIq_subset(const float * const gpu_coor, const float * const gpu_b, const float * const gpu_gv, const float * const gpu_q, float * const gpu_Is, const int Natoms, const int Natoms_subset, const int Nq, const int Ngv, const float scale);
        void wrapper_cudaGV_calcIq_GVVV(const float * const gpu_coor, const float * const gpu_radius, const float * const gpu_gv, const float * const gpu_q, float * const gpu_Is, const int Natoms, const int Nq, const int Ngv, const float sld_solvent);
        void wrapper_cudaGV_calcIq_GVGaussian(const float * const gpu_coor, const float * const gpu_radius, const float * const gpu_gv, const float * const gpu_q, float * const gpu_Is, const int Natoms, const float volume, const int Nq, const int Ngv, const float sld_solvent);
        void wrapper_cudaGV_calcIq_Merge_ND(const float * const gpu_coor, const float * const gpu_b, const float * const gpu_radius, const float * const gpu_mc_coor, const float * const gpu_gv, const float * const gpu_q, float * const gpu_I, const int Natoms, const int N_mc_points, const float volume, const int Nq, const int Ngv, const float sld_solvent);
        void wrapper_cudaGV_calcIq_Merge_NDsubset(const float * const gpu_coor, const float * const gpu_b, const float * const gpu_radius, const float * const gpu_mc_coor, const float * const gpu_gv, const float * const gpu_q, float * const gpu_I, const int Natoms, const int N_mc_points, const int N_mc_points_subset, const float volume, const int Nq, const int Ngv, const float sld_solvent);
        void wrapper_cudaGV_calcIq_Merge_GVVV(const float * const gpu_coor, const float * const gpu_b, const float * const gpu_radius, const float * const gpu_gv, const float * const gpu_q, float * const gpu_Is, const int Natoms, const int Nq, const int Ngv, const float sld_solvent);
        void wrapper_cudaGV_calcIq_Merge_GVGaussian(const float * const gpu_coor, const float * const gpu_b, const float * const gpu_radius, const float * const gpu_gv, const float * const gpu_q, float * const gpu_Is, const int Natoms, const float volume, const int Nq, const int Ngv, const float sld_solvent);

        void wrapper_cudaDebye_calcIq(const float * const gpu_coor, const float * const gpu_b, const float * const gpu_q, float * const gpu_Ia, const int Natoms, const int Nq);
        void wrapper_cudaDebye_calcIq_subset(const float * const gpu_coor, const float * const gpu_b, const float * const gpu_q, float * const gpu_Ia, const int Natoms, const int Natoms_subset, const int Nq, const float scale);
        void wrapper_cudaDebye_calcIq_CustomizedModel(const float * const gpu_coor, const float * const gpu_b, const float * const gpu_radius, const float * const gpu_q, float * const gpu_Ia, const int Natoms, const int Nq);
    }
}

#endif
