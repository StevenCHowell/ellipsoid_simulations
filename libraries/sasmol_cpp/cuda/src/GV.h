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
        int _Ngv;
        const int _Nq;
        const TYPE _qmax;

        TYPE * _gv;
        TYPE * _q;
        TYPE * _Iq;

    // interface
    public:
        inline const TYPE * q() const {return _q;}
        inline const TYPE * Iq() const {return _Iq;}

    // methods
    protected:
        void _calcGV();
        virtual void _setup();

    // meta
    public:
        GV(const int Ngv, const int Nq, const TYPE qmax);
        virtual ~GV();

    // disallow copy and assign
    private:
        inline GV(const GV &);
        inline const GV & operator=(const GV &);
};

#endif
