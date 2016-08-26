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
cudaGV(const int Ngv, const int Nq, const float qmax):
sascalc::GV(Ngv,Nq,qmax)
{
    CALL_CUDA(cudaMalloc((void**)&_gpu_gv, 3*_Ngv*sizeof(float)));
    CALL_CUDA(cudaMalloc((void**)&_gpu_q, _Nq*sizeof(float)));
    CALL_CUDA(cudaMalloc((void**)&_gpu_Is, _Ngv*_Nq*sizeof(float)));
    CALL_CUDA(cudaMallocHost((void**)&_Is, _Ngv*_Nq*sizeof(float)));
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
    CALL_CUDA(cudaMemcpy(_gpu_gv, _gv, _Ngv*3*sizeof(float), cudaMemcpyHostToDevice));
    CALL_CUDA(cudaMemcpy(_gpu_q, _q, _Nq*sizeof(float), cudaMemcpyHostToDevice));
}

//////////////////////////////////////////////////////////
/// calcIq 
//////////////////////////////////////////////////////////
void
sascalc::cudaGV::
calcIq(const float * const gpu_coor, const float * const gpu_b, const int Natoms)
{
    sascalc::cuda::wrapper_cudaGV_calcIq(gpu_coor, gpu_b, _gpu_gv, _gpu_q, _gpu_Is, Natoms, _Nq, _Ngv);
    _harvest();
}

//////////////////////////////////////////////////////////
/// calcIq 
//////////////////////////////////////////////////////////
void
sascalc::cudaGV::
calcIq_subset(const float * const gpu_coor, const float * const gpu_b, const int Natoms, const int Natoms_subset, const float scale)
{
    sascalc::cuda::wrapper_cudaGV_calcIq_subset(gpu_coor, gpu_b, _gpu_gv, _gpu_q, _gpu_Is, Natoms, Natoms_subset, _Nq, _Ngv, scale);
    _harvest();
}

//////////////////////////////////////////////////////////
/// calcIq for the contrast
//////////////////////////////////////////////////////////
void
sascalc::cudaGV::
calcIq_GVVV(const float * const gpu_coor, const float * const gpu_radius, const int Natoms, const float sld_solvent)
{
    sascalc::cuda::wrapper_cudaGV_calcIq_GVVV(gpu_coor, gpu_radius, _gpu_gv, _gpu_q, _gpu_Is, Natoms, _Nq, _Ngv, sld_solvent);
    _harvest();
}

//////////////////////////////////////////////////////////
/// calcIq for the contrast
//////////////////////////////////////////////////////////
void
sascalc::cudaGV::
calcIq_GVGaussian(const float * const gpu_coor, const float * const gpu_radius, const int Natoms, const float volume, const float sld_solvent)
{
    sascalc::cuda::wrapper_cudaGV_calcIq_GVGaussian(gpu_coor, gpu_radius, _gpu_gv, _gpu_q, _gpu_Is, Natoms, volume, _Nq, _Ngv, sld_solvent);
    _harvest();
}

//////////////////////////////////////////////////////////
/// calcIq for the contrast
//////////////////////////////////////////////////////////
void
sascalc::cudaGV::
calcIq_Merge_ND(const float * const gpu_coor, const float * const gpu_b, const float * const gpu_radius, const float * const gpu_mc_coor, const int Natoms, const int N_mc_points, const float volume, const float sld_solvent)
{
    sascalc::cuda::wrapper_cudaGV_calcIq_Merge_ND(gpu_coor, gpu_b, gpu_radius, gpu_mc_coor, _gpu_gv, _gpu_q, _gpu_Is, Natoms, N_mc_points, volume, _Nq, _Ngv, sld_solvent);
    _harvest();
}

//////////////////////////////////////////////////////////
/// calcIq for the contrast
//////////////////////////////////////////////////////////
void
sascalc::cudaGV::
calcIq_Merge_GVVV(const float * const gpu_coor, const float * const gpu_b, const float * const gpu_radius, const int Natoms, const float sld_solvent, const float scale_volume)
{
    sascalc::cuda::wrapper_cudaGV_calcIq_Merge_GVVV(gpu_coor, gpu_b, gpu_radius, _gpu_gv, _gpu_q, _gpu_Is, Natoms, _Nq, _Ngv, sld_solvent, scale_volume);
    _harvest();
}

//////////////////////////////////////////////////////////
/// calcIq for the contrast
//////////////////////////////////////////////////////////
void
sascalc::cudaGV::
calcIq_Merge_GVGaussian(const float * const gpu_coor, const float * const gpu_b, const float * const gpu_radius, const int Natoms, const float volume, const float sld_solvent)
{
    sascalc::cuda::wrapper_cudaGV_calcIq_Merge_GVGaussian(gpu_coor, gpu_b, gpu_radius, _gpu_gv, _gpu_q, _gpu_Is, Natoms, volume, _Nq, _Ngv, sld_solvent);
    _harvest();
}

//////////////////////////////////////////////////////////
/// calcIq for the contrast
//////////////////////////////////////////////////////////
void
sascalc::cudaGV::
calcIq_Merge_NDsubset(const float * const gpu_coor, const float * const gpu_b, const float * const gpu_radius, const float * const gpu_mc_coor, const int Natoms, const int N_mc_points, const int N_mc_points_subset, const float volume, const float sld_solvent)
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
    CALL_CUDA(cudaMemcpy(_Is, _gpu_Is, _Ngv*_Nq*sizeof(float), cudaMemcpyDeviceToHost));
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
