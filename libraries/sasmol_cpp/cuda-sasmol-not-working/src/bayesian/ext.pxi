
# get config
include "config.pxi"

# get common declarations
from altar.common cimport *

# intra-altar references
from altar.problem.ext cimport Problem, XProblem

# abstract base classes
include "AnnealingMethod.pxi"

include "COV.pxi"
include "Metropolis.pxi"
include "SequentialAnnealing.pxi"

IF WITH_MPI:
    include "MPIAnnealing.pxi"

IF WITH_CUDA:
    include "cudaMetropolis.pxi"
    include "CudaAnnealing.pxi"

include "Annealer.pxi"

