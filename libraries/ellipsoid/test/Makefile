LIBS ?= -lellipsoid
INCLUDES ?= -I../include
LIBDIRS ?= -L../lib

INCLUDES += -I/home/hailiang/work/LibrariesSources/

#gsl on MAC
LIBS += -lgsl
LIBS += -lgslcblas
#LIBDIRS += -L/opt/local/lib
#INCLUDES += -I/opt/local/include/gsl
#gsl on entropy
LIBDIRS += -L/share/apps/local/lib
INCLUDES += -I/share/apps/local/include

LDFLAG ?= -std=c++11 -g

test.out: test.cpp ../lib/libellipsoid.a 
	g++ -o $@ $^ $(LDFLAG) $(LIBDIRS) $(LIBS) $(INCLUDES)

.PHONY: clean
clean:
	rm -rf test.out*
