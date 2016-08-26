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
Debye(const int Nq, const TYPE qmax):
_Nq(Nq), _qmax(qmax)
{
    CALL_CUDA(cudaMallocHost((void**)&_q, _Nq*sizeof(TYPE)));
    CALL_CUDA(cudaMallocHost((void**)&_Iq, _Nq*sizeof(TYPE)));
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
calcIq(const TYPE * const coor, const TYPE * const b, const int Natoms)
{
    int i,j,iq;
    TYPE q,I,r;
    TYPE xi,yi,zi,bi,xj,yj,zj,bj;

    for (iq=0; iq<_Nq; ++iq) _Iq[iq] = 0.0;


/*
    // calc I0
    I = 0.0;
    for (i=0; i<Natoms; ++i)
    {
        for (j=i; j<Natoms; ++j) I += b[i]*b[j];
        //std::cout<<i<<" "<<I<<" sofar: "<<_Iq[0]<<std::endl;
        std::cout<<i<<" sofar "<<I<<std::endl;
    }
    _Iq[0] += I;
    //std::cout<<_Iq[0]<<std::endl;
*/

    // calc the rest
    for (iq=0; iq<_Nq; ++iq)
    {
        q = _q[iq];
        for (i=0; i<Natoms; ++i)
        {
            I = 0.0;
            xi = coor[i];
            yi = coor[Natoms + i];
            zi = coor[Natoms*2 + i];
            bi = b[i];
            //I += bi*bi;
            I += bi*bi;
            for (j=i+1; j<Natoms; ++j)
            {
                xj = coor[j];
                yj = coor[Natoms + j];
                zj = coor[Natoms*2 + j];
                bj = b[j];
                r = sqrt(pow((xi-xj),2.0)+pow((yi-yj),2.0)+pow((zi-zj),2.0));
//if (iq==0) std::cout<<"q*r: "<<(q*r)<<" bool: "<<(q*r==0.0)<<" "<<bi<<" "<<bj<<std::endl;
                if (q*r==0.0) I += 2.0*bi*bj;
                else I += 2.0*bi*bj*sin(q*r)/(q*r);
            }
            _Iq[iq] += I;
        }
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
