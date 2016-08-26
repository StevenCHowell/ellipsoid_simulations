//
// Hailiang Zhang
// NIST & UTK
//

#ifndef H_DEBYE
#define H_DEBYE

#include <cuda.h>
#include <cuda_runtime.h>

#include "cudaUtil.h"
#include <sasmol.h>

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
        sasmol::SasMol * _pmol;
        const int _Nq;
        const float _qmax;
        float * _q;
        float * _Iq;

    // interface
    public:
        inline const float * q() const {return _q;}
        inline const float * Iq() const {return _Iq;}

    // methods
    protected:
        virtual void _setup();
    public:
        virtual void calcIq(const float * const gpu_coor, const float * const gpu_b, const int Natoms);

    // meta
    public:
        Debye(const int Nq, const float qmax);
        Debye(sasmol::SasMol * pmol, const int Nq, const float qmax);
        virtual ~Debye();

    // disallow copy and assign
    private:
        inline Debye(const Debye &);
        inline const Debye & operator=(const Debye &);
};

#endif
