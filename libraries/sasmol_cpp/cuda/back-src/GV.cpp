//
// Hailiang Zhang
// NIST & UTK
//

#include "GV.h"
#include <math.h>

//////////////////////////////////////////////////////////
/// constructor
//////////////////////////////////////////////////////////
sascalc::GV::
GV(const int Ngv, const int Nq, const float qmax):
_Ngv(Ngv), _Nq(Nq),_qmax(qmax)
{
    if (_Ngv%2==0)
    {
        std::cout<<"The number of golden vectors should be an odd integer, and so it will be reset to be: "<<_Ngv+1<<std::endl;
        ++_Ngv;
    }
    CALL_CUDA(cudaMallocHost((void**)&_gv, 3*_Ngv*sizeof(float)));
    CALL_CUDA(cudaMallocHost((void**)&_q, _Nq*sizeof(float)));
    CALL_CUDA(cudaMallocHost((void**)&_Iq, _Nq*sizeof(float)));
    _setup();
}

//////////////////////////////////////////////////////////
/// Golden vector calculator
//////////////////////////////////////////////////////////
void
sascalc::GV::
_calcGV()
{
        const float phi_inv = 2.0/(1+sqrt(5.)); // golden ratio
        float cos_theta, sin_theta, phi;
        float qx,qy,qz;
        int igv;
        const int rank = _Ngv/2;
        for (int i=-rank; i<=rank; i++)
        {   
            sin_theta = cos(asin(2.0*i/_Ngv));
            cos_theta = 2.0*i/_Ngv;
            phi = 2*M_PI*i*phi_inv;
            igv = i + rank;
            _gv[igv] = sin_theta*cos(phi);
            _gv[_Ngv+igv] = sin_theta*sin(phi);
            _gv[2*_Ngv+igv] = cos_theta;
        }   
}

//////////////////////////////////////////////////////////
/// set up GV 
//////////////////////////////////////////////////////////
void
sascalc::GV::
_setup()
{
    _calcGV();
    for (int i=0; i<_Nq; ++i) _q[i] = i*(_qmax/(_Nq-1));
}


//////////////////////////////////////////////////////////
/// GV deallocation
//////////////////////////////////////////////////////////
sascalc::GV::
~GV()
{
    CALL_CUDA(cudaFreeHost(_gv));
    CALL_CUDA(cudaFreeHost(_q));
    CALL_CUDA(cudaFreeHost(_Iq));
}
