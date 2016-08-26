# directories
DIRSRC=src/
DIRINC=include/
DIRLIB=lib/
DIRINC_EIGEN=../ # Eigen location

#gsl on MAC
#DIRLIB_GSL=/opt/local/lib
#DIRINC_GSL=/opt/local/include/gsl
#gsl on entropy
#DIRLIB_GSL=/share/apps/local/lib
#DIRINC_GSL=/share/apps/local/include
#gsl on linux
DIRLIB_GSL=/usr/local/lib
DIRINC_GSL=/usr/local/include

# compilers and flags
# must use g++ >= 4.8
#CC=/share/apps/local/bin/g++
CC=g++
CCFLAGS=-fPIC -g -std=c++11

LIBNAME=ellipsoid

LIB=$(DIRLIB)lib$(LIBNAME).a

OBJS=ellipsoid
SOURCES=$(addsuffix .cpp,$(addprefix $(DIRSRC),$(OBJS)))
HEADERS=$(addsuffix .h,$(addprefix $(DIRINC),$(OBJS)))
OBJECTS=$(addsuffix .o,$(addprefix $(DIRLIB),$(OBJS)))

# dependencies
all: $(LIB)

$(LIB): $(OBJECTS)
	ar rv $@ $^

$(OBJECTS): $(DIRLIB)%.o: $(DIRSRC)%.cpp $(DIRINC)%.h
	$(CC) -c $< -o $@ $(CCFLAGS) -I$(DIRINC) -I$(DIRINC_EIGEN) -I$(DIRINC_GSL)


.PHONY: clean
clean:
	rm -rf $(OBJECTS) $(LIB)
