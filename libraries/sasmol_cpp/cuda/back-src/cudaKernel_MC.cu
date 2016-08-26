//
// Hailiang Zhang
// NIST & UTK
//

#include "cudaKernel_MC.h"


/// @par Main functionality
__global__ void setupCurandKernel(curandState * const state, const int seed, const int n)
{
    int id = threadIdx.x + blockIdx.x*blockDim.x;
    if(id<n) curand_init(seed, id, 0, &state[id]);
}

/// @par Main functionality
__global__ void cudaSetupGrid(const float * const gpu_coor, const float * const gpu_radii, int * const gpu_grid, const int Natoms,
                              const float x0, const float y0, const float z0, const float delta, const int Nx, const int Ny, const int Nz)
{
	// get the thread id
    const int idx = blockIdx.x*blockDim.x + threadIdx.x;
    if (idx >= Natoms) return;

    float x,y,z,radius;
    x = gpu_coor[idx];
    y = gpu_coor[Natoms + idx];
    z = gpu_coor[Natoms*2 + idx];
    radius = gpu_radii[idx];

    int ix_low = max(int(((x-radius)-x0)/delta), 0);
    int ix_high = min(int(((x+radius)-x0)/delta), Nx-1);
    int iy_low = max(int(((y-radius)-y0)/delta), 0);
    int iy_high = min(int(((y+radius)-y0)/delta), Ny-1);
    int iz_low = max(int(((z-radius)-z0)/delta), 0);
    int iz_high = min(int(((z+radius)-z0)/delta), Nz-1);

    for (int ix=ix_low; ix<=ix_high; ++ix)
    {
        for (int iy=iy_low; iy<=iy_high; ++iy)
        {
            for (int iz=iz_low; iz<=iz_high; ++iz)
            {
                gpu_grid[iz*(Nx*Ny)+iy*Nx+ix] = 1;
            }
        }
    }
}

/// @par Main functionality
__global__ void cudaMC(const float * const gpu_coor, const float * const gpu_radii, const int * const gpu_grid, float * const gpu_points,
                       const int Natoms, const int Npoints, int * const gpu_Npoints_tried,
                       const float x0, const float y0, const float z0, const float delta, const int Nx, const int Ny, const int Nz,
                       curandState * const state)
{
	// get the thread id
    const int idx = blockIdx.x*blockDim.x + threadIdx.x;
    if (idx >= Npoints) return;

    //if (idx<10) printf("idx: %d,  Natoms: %d Npoints: %d,  x0: %f y0: %f z0: %f,  xl: %f yl: %f zl: %f,  \n",idx,Natoms,Npoints,x0,y0,z0,xl,yl,zl);

    curandState localState = state[idx];

    float x,y,z;
    float r;
    int n_grid;
    int ix,iy,iz;
    const int Nxyz = Nx*Ny*Nz;
    bool flag_found = false;
    while (!flag_found)
    {
        // get a random grid and continue if it is not a candidate grid
        n_grid = int((1.0-curand_uniform(&localState))*(Nxyz)); // curand_uniform will generate random number (0.0,1.0], but I need [0.0,1.0)
        if (gpu_grid[n_grid]==0) continue;

        // update the number of trials
        atomicAdd(gpu_Npoints_tried, 1);

        // get a random point inside the grid
        ix = (n_grid%(Nx*Ny))%Nx;
        iy = (n_grid%(Nx*Ny))/Nx;
        iz = n_grid/(Nx*Ny);
        x = x0 + ix*delta + curand_uniform(&localState)*delta;
        y = y0 + iy*delta + curand_uniform(&localState)*delta;
        z = z0 + iz*delta + curand_uniform(&localState)*delta;

        // check whether this points is inside the molecule or not
        for (int i=0; i<Natoms; ++i)
        {
            r = sqrt(powf((x-gpu_coor[i]),2.0)+powf((y-gpu_coor[Natoms+i]),2.0)+powf((z-gpu_coor[2*Natoms+i]),2.0));
            if (r<gpu_radii[i])
            {
                flag_found = true;
                gpu_points[idx] = x;
                gpu_points[Npoints+idx] = y;
                gpu_points[2*Npoints+idx] = z;
                break;
            }
        }
    }
}
