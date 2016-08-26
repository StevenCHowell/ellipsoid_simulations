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
Debye(const int Nq, const float qmax, sasmol::SasMol * pmol):
_Nq(Nq), _qmax(qmax), _pmol(pmol)
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
    _B['H']=-0.37390;
    _B['D']=0.6671;
    _B['C'] = 0.6646;
    _B['S'] = 0.2847;
    _B['P'] = 0.513;
    _B['N'] = 0.936;
    _B['O'] = 0.5803;
}

//////////////////////////////////////////////////////////
/// calculate Iq with sasmol
//////////////////////////////////////////////////////////
void
sascalc::Debye::
calcIq()
{
    const int natoms = _pmol->_natoms();
    const Eigen::Array<float,Eigen::Dynamic,Eigen::Dynamic> coor = _pmol->_coor();
    
    // get the bs for all the atoms
    float * const b = (float*)malloc(natoms*sizeof(float));
    for (int i=0; i<natoms; ++i) b[i] = _B.at(_pmol->_atom_name()[i][0]);

    calcIq(coor.data(), b, natoms);

    free(b);
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

    for (iq=0; iq<_Nq; ++iq) _Iq[iq] = 0.0;

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
            I += bi*bi;
            for (j=i+1; j<Natoms; ++j)
            {
                xj = coor[j];
                yj = coor[Natoms + j];
                zj = coor[Natoms*2 + j];
                bj = b[j];
                r = sqrt(pow((xi-xj),2.0)+pow((yi-yj),2.0)+pow((zi-zj),2.0));
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
