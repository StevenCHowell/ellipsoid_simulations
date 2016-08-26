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

#include <map>

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
        const float _qmax;
        float * _q;
        float * _Iq;
        sasmol::SasMol * _pmol;
        std::map<char,float> _B;

    // interface
    public:
        inline const float * q() const {return _q;}
        inline const float * Iq() const {return _Iq;}

    // methods
    protected:
        virtual void _setup();
    public:
        virtual void calcIq(const float * const gpu_coor, const float * const gpu_b, const int Natoms);
        virtual void calcIq();

    // meta
    public:
        Debye(const int Nq, const float qmax, sasmol::SasMol * pmol = 0);
        virtual ~Debye();

    // disallow copy and assign
    private:
        inline Debye(const Debye &);
        inline const Debye & operator=(const Debye &);
};

#endif
