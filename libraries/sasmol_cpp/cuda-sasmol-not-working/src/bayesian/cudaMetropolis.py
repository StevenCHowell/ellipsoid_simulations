# -*- coding: utf-8 -*-


from .Metropolis import Metropolis


class cudaMetropolis(Metropolis):

    # private data
    from .ext import XCudaMetropolis as XMetropolis # override my extension type


# end of file 
