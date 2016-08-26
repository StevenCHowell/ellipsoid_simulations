//
// Hailiang Zhang
// NIST & UTK
//

#include "cudaMC.h"
#include "cudaUtil.h"
#include <sasmol.h>
#include <sasio.h>

//////////////////////////////////////////////////////////
/// constructor
//////////////////////////////////////////////////////////
sascalc::cudaMC::
cudaMC(const int Npoints, const int seed):
MC(Npoints,seed)
{
    CALL_CUDA(cudaMalloc((void**)&_gpu_points, 3*_Npoints*sizeof(TYPE)));
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
TYPE
sascalc::cudaMC::
fill(const TYPE * const gpu_coor, const TYPE * const gpu_radii, const int Natoms, const TYPE x0, const TYPE y0, const TYPE z0, const TYPE xl, const TYPE yl, const TYPE zl)
{
    // initialization
    const TYPE delta = 2.0; // ZHl hardwired
    const int Nx = int(xl/delta)+1;
    const int Ny = int(yl/delta)+1;
    const int Nz = int(zl/delta)+1;
    int * gpu_grid;
    CALL_CUDA(cudaMalloc((void**)&gpu_grid, Nx*Ny*Nz*sizeof(int)));
    CALL_CUDA(cudaMemset(gpu_grid, 0, Nx*Ny*Nz*sizeof(int)));
    int * cpu_grid;
    CALL_CUDA(cudaMallocHost((void**)&cpu_grid, Nx*Ny*Nz*sizeof(int)));
    int * gpu_Npoints_tried;
    CALL_CUDA(cudaMalloc((void**)&gpu_Npoints_tried, 1*sizeof(int)));
    CALL_CUDA(cudaMemset(gpu_Npoints_tried, 0, 1*sizeof(int)));
    int * cpu_Npoints_tried;
    CALL_CUDA(cudaMallocHost((void**)&cpu_Npoints_tried, 1*sizeof(int)));

    // calculation
    sascalc::cuda::wrapper_cudaMC_setupGrid(gpu_coor, gpu_radii, gpu_grid, Natoms, x0, y0, z0, delta, Nx, Ny, Nz);
    sascalc::cuda::wrapper_cudaMC_fill(gpu_coor, gpu_radii, gpu_grid, _gpu_points, Natoms, _Npoints, gpu_Npoints_tried, x0,y0,z0,delta,Nx,Ny,Nz,_state);

    // collecting
    _harvest();

    // get volume
    CALL_CUDA(cudaMemcpy(cpu_grid, gpu_grid, Nx*Ny*Nz*sizeof(int), cudaMemcpyDeviceToHost));
    CALL_CUDA(cudaMemcpy(cpu_Npoints_tried, gpu_Npoints_tried, 1*sizeof(int), cudaMemcpyDeviceToHost));
    int count_filled_grid = 0;
    for (int i=0; i<Nx*Ny*Nz; ++i) {if (cpu_grid[i]) ++count_filled_grid;}
    TYPE volume = (delta*delta*delta*count_filled_grid)*(TYPE(_Npoints)/TYPE(*cpu_Npoints_tried));

    // cleaning
    CALL_CUDA(cudaMalloc((void**)&gpu_Npoints_tried, 1*sizeof(int)));
    CALL_CUDA(cudaFree(gpu_grid));
    CALL_CUDA(cudaFree(gpu_Npoints_tried));
    CALL_CUDA(cudaFreeHost(cpu_grid));
    CALL_CUDA(cudaFreeHost(cpu_Npoints_tried));

    // return
    return volume;
}

//////////////////////////////////////////////////////////
/// harvest
//////////////////////////////////////////////////////////
void
sascalc::cudaMC::
_harvest()
{
    CALL_CUDA(cudaMemcpy(_points, _gpu_points, _Npoints*3*sizeof(TYPE), cudaMemcpyDeviceToHost));
}

//////////////////////////////////////////////////////////
/// send 
//////////////////////////////////////////////////////////
void
sascalc::cudaMC::
send()
{
    CALL_CUDA(cudaMemcpy(_gpu_points, _points, _Npoints*3*sizeof(TYPE), cudaMemcpyHostToDevice));
}

//////////////////////////////////////////////////////////
/// destruction
//////////////////////////////////////////////////////////
sascalc::cudaMC::
~cudaMC()
{
    CALL_CUDA(cudaFree(_gpu_points));
}
