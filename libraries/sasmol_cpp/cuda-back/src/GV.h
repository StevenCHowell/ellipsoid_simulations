//
// Hailiang Zhang
// NIST & UTK
//

#ifndef H_GV
#define H_GV

#include <cuda.h>
#include <cuda_runtime.h>

#include "cudaUtil.h"

// declaration
namespace sascalc
{
    /// @brief (Brief description)<br>
    ///
    /// (Detailed descriptions)<br>
    class GV;
}

// definitions
class sascalc::GV
{
    // data
    protected:
        const int _rank;
        const int _NGV;
        const int _Nq;
        const float _qmax;

        float * _gv;
        float * _q;
        float * _Iq;

    // interface
    public:
        inline const float * q() const {return _q;}
        inline const float * Iq() const {return _Iq;}

    // methods
    protected:
        void _calcGV();
        virtual void _setup();

    // meta
    public:
        GV(const int rank, const int Nq, const float qmax);
        virtual ~GV();

    // disallow copy and assign
    private:
        inline GV(const GV &);
        inline const GV & operator=(const GV &);
};

#endif
