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
        float * _gpu_gv;
        float * _gpu_q;
        float * _Ia;
        float * _gpu_Ia;

    // methods
    protected:
        virtual void _setup();
        void _harvest();
    public:
        virtual void calcIq(const float * const gpu_coor, const float * const gpu_b);

    // meta
    public:
        cudaDebye(const int Nq, const float qmax, const int Natoms);
        virtual ~cudaDebye();

    // disallow copy and assign
    private:
        inline cudaDebye(const cudaDebye &);
        inline const cudaDebye & operator=(const cudaDebye &);
};

#endif
