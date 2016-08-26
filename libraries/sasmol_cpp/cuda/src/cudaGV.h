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
        TYPE * _gpu_gv;
        TYPE * _gpu_q;
        TYPE * _Is;
        TYPE * _gpu_Is;

    // methods
    protected:
        virtual void _setup();
        void _harvest();
    public:
        virtual void calcIq(const TYPE * const gpu_coor, const TYPE * const gpu_b, const int Natoms);
        virtual void calcIq_subset(const TYPE * const gpu_coor, const TYPE * const gpu_b, const int Natoms, const int Natoms_subset, const TYPE scale=1.0);
        virtual void calcIq_GVVV(const TYPE * const gpu_coor, const TYPE * const gpu_radius, const int Natoms, const TYPE sld_solvent=-0.005618666666666664);
        virtual void calcIq_GVGaussian(const TYPE * const gpu_coor, const TYPE * const gpu_radius, const int Natoms, const TYPE volume, const TYPE sld_solvent=-0.005618666666666664);
        virtual void calcIq_Merge_ND(const TYPE * const gpu_coor, const TYPE * const gpu_b, const TYPE * const gpu_radius, const TYPE * const gpu_mc_coor, const int Natoms, const int N_mc_points, const TYPE volume, const TYPE sld_solvent=-0.005618666666666664);
        virtual void calcIq_Merge_NDsubset(const TYPE * const gpu_coor, const TYPE * const gpu_b, const TYPE * const gpu_radius, const TYPE * const gpu_mc_coor, const int Natoms, const int N_mc_points, const int N_mc_points_subset, const TYPE volume, const TYPE sld_solvent=-0.005618666666666664);
        virtual void calcIq_Merge_GVVV(const TYPE * const gpu_coor, const TYPE * const gpu_b, const TYPE * const gpu_radius, const int Natoms, const TYPE sld_solvent=-0.005618666666666664, const TYPE scale_volume=0.7);
        virtual void calcIq_Merge_GVGaussian(const TYPE * const gpu_coor, const TYPE * const gpu_b, const TYPE * const gpu_radius, const int Natoms, const TYPE volume, const TYPE sld_solvent=-0.005618666666666664);

    // meta
    public:
        cudaGV(const int Ngv, const int Nq, const TYPE qmax);
        virtual ~cudaGV();

    // disallow copy and assign
    private:
        inline cudaGV(const cudaGV &);
        inline const cudaGV & operator=(const cudaGV &);
};

#endif
