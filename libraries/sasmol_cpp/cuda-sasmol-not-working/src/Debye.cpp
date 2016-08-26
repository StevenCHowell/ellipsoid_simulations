//
// Hailiang Zhang
// NIST & UTK
//

#include "Debye.h"
#include <math.h>

//////////////////////////////////////////////////////////
/// construction
//////////////////////////////////////////////////////////
sascalc::Debye::
Debye(const int Nq, const float qmax):
_Nq(Nq), _qmax(qmax)
{
    CALL_CUDA(cudaMallocHost((void**)&_q, _Nq*sizeof(float)));
    CALL_CUDA(cudaMallocHost((void**)&_Iq, _Nq*sizeof(float)));
    _setup();
}

//////////////////////////////////////////////////////////
/// construction
//////////////////////////////////////////////////////////
sascalc::Debye::
Debye(sasmol::SasMol * const pmol, const int Nq, const float qmax):
_pmol(pmol), _Nq(Nq), _qmax(qmax)
{
    CALL_CUDA(cudaMallocHost((void**)&_q, _Nq*sizeof(float)));
    CALL_CUDA(cudaMallocHost((void**)&_Iq, _Nq*sizeof(float)));
    _setup();
}

//////////////////////////////////////////////////////////
/// set up
//////////////////////////////////////////////////////////
void
sascalc::Debye::
_setup()
{
    for (int i=0; i<_Nq; ++i) _q[i] = i*(_qmax/(_Nq-1));
}

//////////////////////////////////////////////////////////
/// calculate Iq
//////////////////////////////////////////////////////////
void
sascalc::Debye::
calcIq(const float * const coor, const float * const b, const int Natoms)
{
    int i,j,iq;
    float q,I,r;
    float xi,yi,zi,bi,xj,yj,zj,bj;

    // calc I0
    for (i=0; i<Natoms; ++i) for (j=i; j<Natoms; ++j) _Iq[0] += b[i]*b[j];

    // calc the rest
    for (iq=1; iq<_Nq; ++iq)
    {
        q = _q[iq];
        I = 0.0;
        for (i=0; i<Natoms; ++i)
        {
            xi = coor[i];
            yi = coor[Natoms + i];
            zi = coor[Natoms*2 + i];
            bi = b[i];
            I += bi*bi;
            for (j=i+1; j<Natoms; ++j)
            {
                xj = coor[j];
                yj = coor[Natoms + j];
                zj = coor[Natoms*2 + j];
                bj = b[j];
                r = sqrt(pow((xi-xj),2.0)+pow((yi-yj),2.0)+pow((zi-zj),2.0));
                I += bi*bj*sin(q*r)/(q*r);
            }
        }
        _Iq[iq] = I;
    }
}

//////////////////////////////////////////////////////////
/// destruction
//////////////////////////////////////////////////////////
sascalc::Debye::
~Debye()
{
    CALL_CUDA(cudaFreeHost(_q));
    CALL_CUDA(cudaFreeHost(_Iq));
}
