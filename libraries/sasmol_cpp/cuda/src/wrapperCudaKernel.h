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
        void wrapper_cudaMC_setupGrid(const TYPE * const gpu_coor, const TYPE * const gpu_radii, int * const gpu_grid, const int Natoms,
                                      const TYPE x0, const TYPE y0, const TYPE z0, const TYPE delta, const int Nx, const int Ny, const int Nz);
        void wrapper_cudaMC_fill(const TYPE * const gpu_coor, const TYPE * const gpu_radii, const int * const gpu_grid, TYPE * const gpu_points,
                                const int Natoms, const int Npoints, int * const gpu_Npoints_tried,
                                const TYPE x0, const TYPE y0, const TYPE z0, const TYPE delta, const int Nx, const int Ny, const int Nz,
                                curandState * const state);

        void wrapper_cudaGV_calcIq(const TYPE * const gpu_coor, const TYPE * const gpu_b, const TYPE * const gpu_gv, const TYPE * const gpu_q, TYPE * const gpu_Is, const int Natoms, const int Nq, const int Ngv);
        void wrapper_cudaGV_calcIq_subset(const TYPE * const gpu_coor, const TYPE * const gpu_b, const TYPE * const gpu_gv, const TYPE * const gpu_q, TYPE * const gpu_Is, const int Natoms, const int Natoms_subset, const int Nq, const int Ngv, const TYPE scale);
        void wrapper_cudaGV_calcIq_GVVV(const TYPE * const gpu_coor, const TYPE * const gpu_radius, const TYPE * const gpu_gv, const TYPE * const gpu_q, TYPE * const gpu_Is, const int Natoms, const int Nq, const int Ngv, const TYPE sld_solvent);
        void wrapper_cudaGV_calcIq_GVGaussian(const TYPE * const gpu_coor, const TYPE * const gpu_radius, const TYPE * const gpu_gv, const TYPE * const gpu_q, TYPE * const gpu_Is, const int Natoms, const TYPE volume, const int Nq, const int Ngv, const TYPE sld_solvent);
        void wrapper_cudaGV_calcIq_Merge_ND(const TYPE * const gpu_coor, const TYPE * const gpu_b, const TYPE * const gpu_radius, const TYPE * const gpu_mc_coor, const TYPE * const gpu_gv, const TYPE * const gpu_q, TYPE * const gpu_I, const int Natoms, const int N_mc_points, const TYPE volume, const int Nq, const int Ngv, const TYPE sld_solvent);
        void wrapper_cudaGV_calcIq_Merge_NDsubset(const TYPE * const gpu_coor, const TYPE * const gpu_b, const TYPE * const gpu_radius, const TYPE * const gpu_mc_coor, const TYPE * const gpu_gv, const TYPE * const gpu_q, TYPE * const gpu_I, const int Natoms, const int N_mc_points, const int N_mc_points_subset, const TYPE volume, const int Nq, const int Ngv, const TYPE sld_solvent);
        void wrapper_cudaGV_calcIq_Merge_GVVV(const TYPE * const gpu_coor, const TYPE * const gpu_b, const TYPE * const gpu_radius, const TYPE * const gpu_gv, const TYPE * const gpu_q, TYPE * const gpu_Is, const int Natoms, const int Nq, const int Ngv, const TYPE sld_solvent, const TYPE scale_volume);
        void wrapper_cudaGV_calcIq_Merge_GVGaussian(const TYPE * const gpu_coor, const TYPE * const gpu_b, const TYPE * const gpu_radius, const TYPE * const gpu_gv, const TYPE * const gpu_q, TYPE * const gpu_Is, const int Natoms, const TYPE volume, const int Nq, const int Ngv, const TYPE sld_solvent);

        void wrapper_cudaDebye_calcIq(const TYPE * const gpu_coor, const TYPE * const gpu_b, const TYPE * const gpu_q, TYPE * const gpu_Ia, const int Natoms, const int Nq);
        void wrapper_cudaDebye_calcIq_subset(const TYPE * const gpu_coor, const TYPE * const gpu_b, const TYPE * const gpu_q, TYPE * const gpu_Ia, const int Natoms, const int Natoms_subset, const int Nq, const TYPE scale);
        void wrapper_cudaDebye_calcIq_CustomizedModel(const TYPE * const gpu_coor, const TYPE * const gpu_b, const TYPE * const gpu_radius, const TYPE * const gpu_q, TYPE * const gpu_Ia, const int Natoms, const int Nq);
    }
}

#endif
