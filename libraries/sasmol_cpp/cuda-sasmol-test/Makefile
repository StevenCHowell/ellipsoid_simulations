# directories
DIRSRC=src/
DIRLIB=lib/

INCDIR=-I/home/hailiang/work/sascalc/sasmol_cpp/include -I/usr/local/cuda-6.0/include

# compilers and flags
CC=g++
CCFLAGS=-fPIC -g

NVCC=/usr/local/cuda-6.0/bin/nvcc
NVCCFLAGS=-Xcompiler -fPIC -DCUCC_ID="nVidia" -DBLD_PROCEDURE -DWITH_GSL -DHAVE_INLINE -DWITH_CUDA -DWITH_MPI -DPARALLEL -DWITH_PARALLEL_IO  $(CUDA_ARCH) -g -G -arch=sm_30

LIBNAME=sascalc
CUDA_LIBNAME=cudaSascalc
EXTLIBNAME=sasmol

LIB=$(DIRLIB)lib$(LIBNAME).a
CUDA_LIB=$(DIRLIB)lib$(CUDA_LIBNAME).a

OBJS=GV MC Debye cudaMC cudaGV cudaDebye
SOURCES=$(addsuffix .cpp,$(addprefix $(DIRSRC),$(OBJS)))
HEADERS=$(SOURCES:.cpp=.h)
OBJECTS=$(addsuffix .o,$(addprefix $(DIRLIB),$(OBJS)))

CUDA_OBJS=cudaKernel_MC cudaKernel_GV cudaKernel_Debye wrapperCudaKernel
CUDA_SOURCES=$(addsuffix .cu,$(addprefix $(DIRSRC),$(CUDA_OBJS)))
CUDA_HEADERS=$(CUDA_SOURCES:.cu=.h)
CUDA_OBJECTS=$(addsuffix .o,$(addprefix $(DIRLIB),$(CUDA_OBJS)))

# dependencies
all: $(LIB) $(CUDA_LIB)

$(LIB): $(OBJECTS)
	ar rv $@ $^

$(OBJECTS): $(DIRLIB)%.o: $(DIRSRC)%.cpp $(DIRSRC)%.h
	$(CC) -c $< -o $@ $(CCFLAGS) $(INCDIR)

$(CUDA_LIB): $(CUDA_OBJECTS)
	ar rv $@ $^

$(CUDA_OBJECTS): $(DIRLIB)%.o: $(DIRSRC)%.cu $(DIRSRC)%.h
	$(NVCC) -c $< -o $@ $(NVCCFLAGS)


.PHONY: clean
clean:
	rm -rf $(OBJECTS) $(CUDA_OBJECTS) $(LIB) $(CUDA_LIB)
