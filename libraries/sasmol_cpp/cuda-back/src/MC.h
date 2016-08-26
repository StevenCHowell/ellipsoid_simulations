//
// Hailiang Zhang
// NIST & UTK
//

#ifndef H_MC
#define H_MC

#include <cuda.h>
#include <cuda_runtime.h>
#include "cudaUtil.h"

// declaration
namespace sascalc
{
    /// @brief (Brief description)<br>
    ///
    /// (Detailed descriptions)<br>
    class MC;
}

// definitions
class sascalc::MC
{
    // data
    protected:
        const int _Npoints; // the number of monter-carlo points to be generated
        float * _points;
        int _seed;

    // interface
    public:
        inline int Npoints() const {return _Npoints;}
        inline const float * points() const {return _points;}

    // methods
    protected:
        virtual void fill(){}// not implemented yet

    // meta
    public:
        MC(const int Npoints, const int seed);
        virtual ~MC();

    // disallow copy and assign
    private:
        inline MC(const MC &);
        inline const MC & operator=(const MC &);
};

#endif
