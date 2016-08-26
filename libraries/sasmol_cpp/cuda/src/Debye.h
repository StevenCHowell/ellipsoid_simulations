//
// Hailiang Zhang
// NIST & UTK
//

#ifndef H_DEBYE
#define H_DEBYE

#include <cuda.h>
#include <cuda_runtime.h>

#include "cudaUtil.h"

// declaration
namespace sascalc
{
    /// @brief (Brief description)<br>
    ///
    /// (Detailed descriptions)<br>
    class Debye;
}

// definitions
class sascalc::Debye
{
    // data
    protected:
        const int _Nq;
        const TYPE _qmax;
        TYPE * _q;
        TYPE * _Iq;

    // interface
    public:
        inline const TYPE * q() const {return _q;}
        inline const TYPE * Iq() const {return _Iq;}

    // methods
    protected:
        virtual void _setup();
    public:
        virtual void calcIq(const TYPE * const gpu_coor, const TYPE * const gpu_b, const int Natoms);

    // meta
    public:
        Debye(const int Nq, const TYPE qmax);
        virtual ~Debye();

    // disallow copy and assign
    private:
        inline Debye(const Debye &);
        inline const Debye & operator=(const Debye &);
};

#endif
