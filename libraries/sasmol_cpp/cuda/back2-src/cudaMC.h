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
        float * _gpu_points;
        curandState * _state;

    // interface
    public:
        inline const float * gpu_points() const {return _gpu_points;}

    // methods
    protected:
        void _setupRand();
        void _harvest();
    public:
        float fill(const float * const gpu_coor, const float * const gpu_radii, const int Natoms, const float x0, const float y0, const float z0, const float xl, const float yl, const float zl);
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
