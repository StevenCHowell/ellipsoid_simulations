//
// Hailiang Zhang
// NIST & UTK
//

#include "cudaMC.h"

//////////////////////////////////////////////////////////
/// constructor
//////////////////////////////////////////////////////////
sascalc::cudaMC::
cudaMC(const int Npoints, const int seed):
MC(Npoints,seed)
{
    CALL_CUDA(cudaMalloc((void**)&_gpu_points, 3*_Npoints*sizeof(float)));
    _setupRand();
}

//////////////////////////////////////////////////////////
/// set up the random generator
//////////////////////////////////////////////////////////
void
sascalc::cudaMC::
_setupRand()
{
    CALL_CUDA(cudaMalloc((void **)&_state, _Npoints*sizeof(curandState)));
    sascalc::cuda::wrapper_cudaMC_setupCurandState(_state, _seed, _Npoints);
}

//////////////////////////////////////////////////////////
/// fill the monte carlo points
//////////////////////////////////////////////////////////
void
sascalc::cudaMC::
fill(const float * const gpu_coor, const float * const gpu_radii, const int Natoms, const float x0, const float y0, const float z0, const float xl, const float yl, const float zl)
{
    sascalc::cuda::wrapper_cudaMC_fill(gpu_coor, gpu_radii, _gpu_points, Natoms, _Npoints,x0,y0,z0,xl,yl,zl,_state);
    _harvest();
}

//////////////////////////////////////////////////////////
/// harvest
//////////////////////////////////////////////////////////
void
sascalc::cudaMC::
_harvest()
{
    CALL_CUDA(cudaMemcpy(_points, _gpu_points, _Npoints*3*sizeof(float), cudaMemcpyDeviceToHost));
}

//////////////////////////////////////////////////////////
/// destruction
//////////////////////////////////////////////////////////
sascalc::cudaMC::
~cudaMC()
{
    CALL_CUDA(cudaFree(_gpu_points));
}
