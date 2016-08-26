//
// Hailiang Zhang
// NIST & UTK
//

#ifndef H_CUDAMC
#define H_CUDAMC

#include "MC.h"
#include "wrapperCudaKernel.h"

// declaration
namespace sascalc
{
    /// @brief (Brief description)<br>
    ///
    /// (Detailed descriptions)<br>
    class cudaMC;
}

// definitions
class sascalc::cudaMC:
    public sascalc::MC
{
    // data
    protected:
        TYPE * _gpu_points;
        curandState * _state;

    // interface
    public:
        inline const TYPE * gpu_points() const {return _gpu_points;}

    // methods
    protected:
        void _setupRand();
        void _harvest();
    public:
        TYPE fill(const TYPE * const gpu_coor, const TYPE * const gpu_radii, const int Natoms, const TYPE x0, const TYPE y0, const TYPE z0, const TYPE xl, const TYPE yl, const TYPE zl);
        void send();

    // meta
    public:
        cudaMC(const int Npoints, const int seed);
        virtual ~cudaMC();

    // disallow copy and assign
    private:
        inline cudaMC(const cudaMC &);
        inline const cudaMC & operator=(const cudaMC &);
};

#endif
