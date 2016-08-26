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
	sascalc::cuda::wrapper_cudaDebye_calcIq(gpu_coor, gpu_b, _gpu_q, _gpu_Ia, _Natoms, _Nq);
    _harvest();
}

//////////////////////////////////////////////////////////
/// calcIq for a customized model
//////////////////////////////////////////////////////////
void
sascalc::cudaDebye::
calcIq_CustomizedModel(const float * const gpu_coor, const float * const gpu_b, const float * const gpu_radius)
{
	sascalc::cuda::wrapper_cudaDebye_calcIq_CustomizedModel(gpu_coor, gpu_b, gpu_radius, _gpu_q, _gpu_Ia, _Natoms, _Nq);
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
