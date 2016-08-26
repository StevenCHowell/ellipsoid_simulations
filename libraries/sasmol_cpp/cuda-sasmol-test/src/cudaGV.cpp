//
// Hailiang Zhang
// NIST & UTK
//

#include "cudaGV.h"
#include <math.h>

//////////////////////////////////////////////////////////
/// construction
//////////////////////////////////////////////////////////
sascalc::cudaGV::
cudaGV(const int rank, const int Nq, const float qmax):
sascalc::GV(rank,Nq,qmax)
{
    CALL_CUDA(cudaMalloc((void**)&_gpu_gv, 3*_NGV*sizeof(float)));
    CALL_CUDA(cudaMalloc((void**)&_gpu_q, _Nq*sizeof(float)));
    CALL_CUDA(cudaMalloc((void**)&_gpu_Is, _NGV*_Nq*sizeof(float)));
    CALL_CUDA(cudaMallocHost((void**)&_Is, _NGV*_Nq*sizeof(float)));
    _setup();
}

//////////////////////////////////////////////////////////
/// set up 
//////////////////////////////////////////////////////////
void
sascalc::cudaGV::
_setup()
{
    GV::_setup();
    CALL_CUDA(cudaMemcpy(_gpu_gv, _gv, _NGV*3*sizeof(float), cudaMemcpyHostToDevice));
    CALL_CUDA(cudaMemcpy(_gpu_q, _q, _Nq*sizeof(float), cudaMemcpyHostToDevice));
}

//////////////////////////////////////////////////////////
/// calcIq 
//////////////////////////////////////////////////////////
void
sascalc::cudaGV::
calcIq(const float * const gpu_coor, const float * const gpu_b, const int Natoms)
{
    sascalc::cuda::wrapper_cudaGV_calcIq(gpu_coor, gpu_b, _gpu_gv, _gpu_q, _gpu_Is, Natoms, _Nq, _NGV);
    _harvest();
}

//////////////////////////////////////////////////////////
/// harvest
//////////////////////////////////////////////////////////
void
sascalc::cudaGV::
_harvest()
{
    CALL_CUDA(cudaMemcpy(_Is, _gpu_Is, _NGV*_Nq*sizeof(float), cudaMemcpyDeviceToHost));
    for (int i=0; i<_Nq; ++i)
    {
        _Iq[i] = 0.0;
        for (int j=0; j<_NGV; ++j)
        {
            _Iq[i] += _Is[j*_Nq+i];
        }
        _Iq[i] /= _NGV;
    }
}

//////////////////////////////////////////////////////////
/// destruction
//////////////////////////////////////////////////////////
sascalc::cudaGV::
~cudaGV()
{
    CALL_CUDA(cudaFree(_gpu_gv));
    CALL_CUDA(cudaFree(_gpu_q));
    CALL_CUDA(cudaFree(_gpu_Is));
    CALL_CUDA(cudaFreeHost(_Is));
}
