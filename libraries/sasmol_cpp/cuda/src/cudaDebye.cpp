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
cudaDebye(const int Nq, const TYPE qmax, const int Natoms):
Debye(Nq,qmax),
_Natoms(Natoms)
{
    CALL_CUDA(cudaMallocHost((void**)&_Ia, Natoms*Nq*sizeof(TYPE)));
    CALL_CUDA(cudaMalloc((void**)&_gpu_q, Nq*sizeof(TYPE)));
    CALL_CUDA(cudaMalloc((void**)&_gpu_Ia, Natoms*Nq*sizeof(TYPE)));
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
    CALL_CUDA(cudaMemcpy(_gpu_q, _q, _Nq*sizeof(TYPE), cudaMemcpyHostToDevice));
}

//////////////////////////////////////////////////////////
/// calcIq 
//////////////////////////////////////////////////////////
void
sascalc::cudaDebye::
calcIq(const TYPE * const gpu_coor, const TYPE * const gpu_b)
{
	sascalc::cuda::wrapper_cudaDebye_calcIq(gpu_coor, gpu_b, _gpu_q, _gpu_Ia, _Natoms, _Nq);
    _harvest();
}

//////////////////////////////////////////////////////////
/// calcIq 
//////////////////////////////////////////////////////////
void
sascalc::cudaDebye::
calcIq_subset(const TYPE * const gpu_coor, const TYPE * const gpu_b, const int Natoms_subset, const TYPE scale)
{
	sascalc::cuda::wrapper_cudaDebye_calcIq_subset(gpu_coor, gpu_b, _gpu_q, _gpu_Ia, _Natoms, Natoms_subset, _Nq, scale);
    _harvest_subset(Natoms_subset);
}

//////////////////////////////////////////////////////////
/// calcIq for a customized model
//////////////////////////////////////////////////////////
void
sascalc::cudaDebye::
calcIq_CustomizedModel(const TYPE * const gpu_coor, const TYPE * const gpu_b, const TYPE * const gpu_radius)
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
    CALL_CUDA(cudaMemcpy(_Ia, _gpu_Ia, _Natoms*_Nq*sizeof(TYPE), cudaMemcpyDeviceToHost));
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
/// harvest Debye for subset
//////////////////////////////////////////////////////////
void
sascalc::cudaDebye::
_harvest_subset(const int Natoms_subset)
{
    CALL_CUDA(cudaMemcpy(_Ia, _gpu_Ia, _Natoms*_Nq*sizeof(TYPE), cudaMemcpyDeviceToHost));
    for (int i=0; i<_Nq; ++i)
    {
        _Iq[i] = 0.0;
        for (int j=0; j<Natoms_subset; ++j)
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
