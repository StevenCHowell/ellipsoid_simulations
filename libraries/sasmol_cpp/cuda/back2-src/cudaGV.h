//
// Hailiang Zhang
// NIST & UTK
//

#ifndef H_CUDA_CUDAGV
#define H_CUDA_CUDAGV

#include <cuda.h>
#include <cuda_runtime.h>

#include "GV.h"
#include "cudaUtil.h"
#include "wrapperCudaKernel.h"

// declaration
namespace sascalc
{
    /// @brief (Brief description)<br>
    ///
    /// (Detailed descriptions)<br>
    class cudaGV;
}

// definitions
class sascalc::cudaGV :
    public sascalc::GV
{
    // data
    protected:
        float * _gpu_gv;
        float * _gpu_q;
        float * _Is;
        float * _gpu_Is;

    // methods
    protected:
        virtual void _setup();
        void _harvest();
    public:
        virtual void calcIq(const float * const gpu_coor, const float * const gpu_b, const int Natoms);
        virtual void calcIq_subset(const float * const gpu_coor, const float * const gpu_b, const int Natoms, const int Natoms_subset, const float scale=1.0);
        virtual void calcIq_GVVV(const float * const gpu_coor, const float * const gpu_radius, const int Natoms, const float sld_solvent=-0.005618666666666664);
        virtual void calcIq_GVGaussian(const float * const gpu_coor, const float * const gpu_radius, const int Natoms, const float volume, const float sld_solvent=-0.005618666666666664);
        virtual void calcIq_Merge_ND(const float * const gpu_coor, const float * const gpu_b, const float * const gpu_radius, const float * const gpu_mc_coor, const int Natoms, const int N_mc_points, const float volume, const float sld_solvent=-0.005618666666666664);
        virtual void calcIq_Merge_NDsubset(const float * const gpu_coor, const float * const gpu_b, const float * const gpu_radius, const float * const gpu_mc_coor, const int Natoms, const int N_mc_points, const int N_mc_points_subset, const float volume, const float sld_solvent=-0.005618666666666664);
        virtual void calcIq_Merge_GVVV(const float * const gpu_coor, const float * const gpu_b, const float * const gpu_radius, const int Natoms, const float sld_solvent=-0.005618666666666664);
        virtual void calcIq_Merge_GVGaussian(const float * const gpu_coor, const float * const gpu_b, const float * const gpu_radius, const int Natoms, const float volume, const float sld_solvent=-0.005618666666666664);

    // meta
    public:
        cudaGV(const int Ngv, const int Nq, const float qmax);
        virtual ~cudaGV();

    // disallow copy and assign
    private:
        inline cudaGV(const cudaGV &);
        inline const cudaGV & operator=(const cudaGV &);
};

#endif
