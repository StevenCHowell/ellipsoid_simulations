//
// Hailiang Zhang
// NIST & UTK
//

#include "cudaDebye.h"
#include <math.h>

//////////////////////////////////////////////////////////
/// construction
//////////////////////////////////////////////////////////
sascalc::cudaDebye::
cudaDebye(const int Nq, const float qmax, const int Natoms):
Debye(Nq,qmax),
_Natoms(Natoms)
{
    CALL_CUDA(cudaMallocHost((void**)&_Ia, Natoms*Nq*sizeof(float)));
    CALL_CUDA(cudaMalloc((void**)&_gpu_q, Nq*sizeof(float)));
    CALL_CUDA(cudaMalloc((void**)&_gpu_Ia, Natoms*Nq*sizeof(float)));
    _setup();
}

//////////////////////////////////////////////////////////
/// construction
//////////////////////////////////////////////////////////
sascalc::cudaDebye::
cudaDebye(sasmol::SasMol * pmol, const int Nq, const float qmax):
Debye(pmol, Nq,qmax),
_Natoms(pmol->_natoms())
{
    CALL_CUDA(cudaMallocHost((void**)&_Ia, _Natoms*_Nq*sizeof(float)));
    CALL_CUDA(cudaMalloc((void**)&_gpu_q, _Nq*sizeof(float)));
    CALL_CUDA(cudaMalloc((void**)&_gpu_Ia, _Natoms*_Nq*sizeof(float)));
    _setup();
}

//////////////////////////////////////////////////////////
/// set up
//////////////////////////////////////////////////////////
void
sascalc::cudaDebye::
_setup()
{
    Debye::_setup();
    CALL_CUDA(cudaMemcpy(_gpu_q, _q, _Nq*sizeof(float), cudaMemcpyHostToDevice));
}

//////////////////////////////////////////////////////////
/// calcIq 
//////////////////////////////////////////////////////////
void
sascalc::cudaDebye::
calcIq(const float * const gpu_coor, const float * const gpu_b)
{
    wrapperCuda_calcIq(gpu_coor, gpu_b);
    _harvest();
}

//////////////////////////////////////////////////////////
/// harvest Debye
//////////////////////////////////////////////////////////
void
sascalc::cudaDebye::
_harvest()
{
    CALL_CUDA(cudaMemcpy(_Ia, _gpu_Ia, _Natoms*_Nq*sizeof(float), cudaMemcpyDeviceToHost));
    for (int i=0; i<_Nq; ++i)
    {
        _Iq[i] = 0.0;
        for (int j=0; j<_Natoms; ++j)
        {
            _Iq[i] += _Ia[i*_Natoms+j];
        }
    }
}

//////////////////////////////////////////////////////////
/// destruction
//////////////////////////////////////////////////////////
sascalc::cudaDebye::
~cudaDebye()
{
    CALL_CUDA(cudaFreeHost(_Ia));
    CALL_CUDA(cudaFree(_gpu_q));
    CALL_CUDA(cudaFree(_gpu_Ia));
}
