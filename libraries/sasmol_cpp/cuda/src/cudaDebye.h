//
// Hailiang Zhang
// NIST & UTK
//

#ifndef H_CUDA_CUDADEBYE
#define H_CUDA_CUDADEBYE

#include <cuda.h>
#include <cuda_runtime.h>

#include "Debye.h"
#include "cudaUtil.h"
#include "wrapperCudaKernel.h"

// declaration
namespace sascalc
{
    /// @brief (Brief description)<br>
    ///
    /// (Detailed descriptions)<br>
    class cudaDebye;
}

// definitions
class sascalc::cudaDebye :
    public sascalc::Debye
{
    // data
    protected:
        const int _Natoms;
        TYPE * _gpu_gv;
        TYPE * _gpu_q;
        TYPE * _Ia;
        TYPE * _gpu_Ia;

    // methods
    protected:
        virtual void _setup();
        void _harvest();
        void _harvest_subset(const int Natoms_subset);
    public:
        virtual void calcIq(const TYPE * const gpu_coor, const TYPE * const gpu_b);
        virtual void calcIq_subset(const TYPE * const gpu_coor, const TYPE * const gpu_b, const int Natoms_subset, const TYPE scale=1.0);
        virtual void calcIq_CustomizedModel(const TYPE * const gpu_coor, const TYPE * const gpu_b, const TYPE * const gpu_radius);

    // meta
    public:
        cudaDebye(const int Nq, const TYPE qmax, const int Natoms);
        virtual ~cudaDebye();

    // disallow copy and assign
    private:
        inline cudaDebye(const cudaDebye &);
        inline const cudaDebye & operator=(const cudaDebye &);
};

#endif
