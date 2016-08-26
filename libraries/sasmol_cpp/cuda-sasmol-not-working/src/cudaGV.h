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
        void wrapperCuda_calcIq(const float * const gpu_coor, const float * const gpu_b, const int Natoms);

    // meta
    public:
        cudaGV(const int rank, const int Nq, const float qmax);
        cudaGV(sasmol::SasMol * pmol, const int rank, const int Nq, const float qmax);
        virtual ~cudaGV();

    // disallow copy and assign
    private:
        inline cudaGV(const cudaGV &);
        inline const cudaGV & operator=(const cudaGV &);
};

#endif
