0. General
The present package includes the libraries and MC simulation packages for various CG models:
.
|-- libraries
|   |-- ellipsoid
|   `-- sasmol_cpp
`-- simulations
    |-- hard_sphere_monomer
    |-- hard_triellipsoid_flexible
    |-- hard_trisphere_rigid
    |-- LJ_sphere_monomer
    `-- soft_triellipsoid_flexible

Please following the first steps to build the libraries, simulation packages, and run the simulation

1. Build libraries by typing:

"""
cd libraries/sasmol_cpp/
make
cd ../../

cd libraries/ellipsoid
make
cd ../../
"""

NOTE: The location of boost/gsl headers/libraries are hardwired for the entropy head-node. Compilation failure on different machine is likely due to different boost/gsl locations, and modifying the corresponding folder in the Makefile generally fixes this problem.

2. Go to the simulation folder whose name indicates the CG model types (e.g. "soft_triellipsoid_flexible" for soft tri-ellipsoid with flexible linker), modify the "nptmc.h" file to adjust the simulation conditions, and "ellipsoid.h" where applicable for ellipsoid parameters.

3. Compile the simulation program by typing:

"""
make
"""

4. Run the simulation by typing:

"""
./a.out
"""
