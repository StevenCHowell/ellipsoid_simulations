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
cudaGV(const int Ngv, const int Nq, const TYPE qmax):
sascalc::GV(Ngv,Nq,qmax)
{
    CALL_CUDA(cudaMalloc((void**)&_gpu_gv, 3*_Ngv*sizeof(TYPE)));
    CALL_CUDA(cudaMalloc((void**)&_gpu_q, _Nq*sizeof(TYPE)));
    CALL_CUDA(cudaMalloc((void**)&_gpu_Is, _Ngv*_Nq*sizeof(TYPE)));
    CALL_CUDA(cudaMallocHost((void**)&_Is, _Ngv*_Nq*sizeof(TYPE)));
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
    CALL_CUDA(cudaMemcpy(_gpu_gv, _gv, _Ngv*3*sizeof(TYPE), cudaMemcpyHostToDevice));
    CALL_CUDA(cudaMemcpy(_gpu_q, _q, _Nq*sizeof(TYPE), cudaMemcpyHostToDevice));
}

//////////////////////////////////////////////////////////
/// calcIq 
//////////////////////////////////////////////////////////
void
sascalc::cudaGV::
calcIq(const TYPE * const gpu_coor, const TYPE * const gpu_b, const int Natoms)
{
    sascalc::cuda::wrapper_cudaGV_calcIq(gpu_coor, gpu_b, _gpu_gv, _gpu_q, _gpu_Is, Natoms, _Nq, _Ngv);
    _harvest();
}

//////////////////////////////////////////////////////////
/// calcIq 
//////////////////////////////////////////////////////////
void
sascalc::cudaGV::
calcIq_subset(const TYPE * const gpu_coor, const TYPE * const gpu_b, const int Natoms, const int Natoms_subset, const TYPE scale)
{
    sascalc::cuda::wrapper_cudaGV_calcIq_subset(gpu_coor, gpu_b, _gpu_gv, _gpu_q, _gpu_Is, Natoms, Natoms_subset, _Nq, _Ngv, scale);
    _harvest();
}

//////////////////////////////////////////////////////////
/// calcIq for the contrast
//////////////////////////////////////////////////////////
void
sascalc::cudaGV::
calcIq_GVVV(const TYPE * const gpu_coor, const TYPE * const gpu_radius, const int Natoms, const TYPE sld_solvent)
{
    sascalc::cuda::wrapper_cudaGV_calcIq_GVVV(gpu_coor, gpu_radius, _gpu_gv, _gpu_q, _gpu_Is, Natoms, _Nq, _Ngv, sld_solvent);
    _harvest();
}

//////////////////////////////////////////////////////////
/// calcIq for the contrast
//////////////////////////////////////////////////////////
void
sascalc::cudaGV::
calcIq_GVGaussian(const TYPE * const gpu_coor, const TYPE * const gpu_radius, const int Natoms, const TYPE volume, const TYPE sld_solvent)
{
    sascalc::cuda::wrapper_cudaGV_calcIq_GVGaussian(gpu_coor, gpu_radius, _gpu_gv, _gpu_q, _gpu_Is, Natoms, volume, _Nq, _Ngv, sld_solvent);
    _harvest();
}

//////////////////////////////////////////////////////////
/// calcIq for the contrast
//////////////////////////////////////////////////////////
void
sascalc::cudaGV::
calcIq_Merge_ND(const TYPE * const gpu_coor, const TYPE * const gpu_b, const TYPE * const gpu_radius, const TYPE * const gpu_mc_coor, const int Natoms, const int N_mc_points, const TYPE volume, const TYPE sld_solvent)
{
    sascalc::cuda::wrapper_cudaGV_calcIq_Merge_ND(gpu_coor, gpu_b, gpu_radius, gpu_mc_coor, _gpu_gv, _gpu_q, _gpu_Is, Natoms, N_mc_points, volume, _Nq, _Ngv, sld_solvent);
    _harvest();
}

//////////////////////////////////////////////////////////
/// calcIq for the contrast
//////////////////////////////////////////////////////////
void
sascalc::cudaGV::
calcIq_Merge_GVVV(const TYPE * const gpu_coor, const TYPE * const gpu_b, const TYPE * const gpu_radius, const int Natoms, const TYPE sld_solvent, const TYPE scale_volume)
{
    sascalc::cuda::wrapper_cudaGV_calcIq_Merge_GVVV(gpu_coor, gpu_b, gpu_radius, _gpu_gv, _gpu_q, _gpu_Is, Natoms, _Nq, _Ngv, sld_solvent, scale_volume);
    _harvest();
}

//////////////////////////////////////////////////////////
/// calcIq for the contrast
//////////////////////////////////////////////////////////
void
sascalc::cudaGV::
calcIq_Merge_GVGaussian(const TYPE * const gpu_coor, const TYPE * const gpu_b, const TYPE * const gpu_radius, const int Natoms, const TYPE volume, const TYPE sld_solvent)
{
    sascalc::cuda::wrapper_cudaGV_calcIq_Merge_GVGaussian(gpu_coor, gpu_b, gpu_radius, _gpu_gv, _gpu_q, _gpu_Is, Natoms, volume, _Nq, _Ngv, sld_solvent);
    _harvest();
}

//////////////////////////////////////////////////////////
/// calcIq for the contrast
//////////////////////////////////////////////////////////
void
sascalc::cudaGV::
calcIq_Merge_NDsubset(const TYPE * const gpu_coor, const TYPE * const gpu_b, const TYPE * const gpu_radius, const TYPE * const gpu_mc_coor, const int Natoms, const int N_mc_points, const int N_mc_points_subset, const TYPE volume, const TYPE sld_solvent)
{
    sascalc::cuda::wrapper_cudaGV_calcIq_Merge_NDsubset(gpu_coor, gpu_b, gpu_radius, gpu_mc_coor, _gpu_gv, _gpu_q, _gpu_Is, Natoms, N_mc_points, N_mc_points_subset, volume, _Nq, _Ngv, sld_solvent);
    _harvest();
}

//////////////////////////////////////////////////////////
/// harvest
//////////////////////////////////////////////////////////
void
sascalc::cudaGV::
_harvest()
{
    CALL_CUDA(cudaMemcpy(_Is, _gpu_Is, _Ngv*_Nq*sizeof(TYPE), cudaMemcpyDeviceToHost));
    for (int i=0; i<_Nq; ++i)
    {
        _Iq[i] = 0.0;
        for (int j=0; j<_Ngv; ++j)
        {
            _Iq[i] += _Is[j*_Nq+i];
        }
        _Iq[i] /= _Ngv;
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
