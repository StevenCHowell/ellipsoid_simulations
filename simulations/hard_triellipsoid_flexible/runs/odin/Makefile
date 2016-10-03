# created by Hailiang Zhang
# modified by Steven C Howell
# 24 August 2016

EXE=hte.bin

# DIRINC=../../libraries/sasmol_cpp/include/
# DIRLIB=../../libraries/sasmol_cpp/lib/
DIRINC=/usr/local/sasmol_cpp/include/
DIRLIB=/usr/local/sasmol_cpp/lib/

LIB=sasmol

#boost on MAC
#DIRLIB_BOOST=/usr/local/lib/

#boost on Entropy
#DIRLIB_BOOST=/share/apps/local/boost_1_55_0/libs/

#boost on linux
DIRINC_BOOST=/usr/local/lib/
DIRLIB_BOOST=/usr/local/lib/

#CC=/share/apps/local/bin/g++
CC=g++
CCFLAGS=-std=c++11 -c -fPIC -w -g -I$(DIRINC) -I.
LDFLAGS=-std=c++11 -g -w -I$(DIRINC) -I.

LIBRARY=$(DIRLIB)lib$(LIB).a

$(EXE): ../../npt_mc.cpp
	$(CC) $^ $(LDFLAGS) -L$(DIRLIB) -l$(LIB) -L$(DIRLIB) -L$(DIRLIB_BOOST) -I$(DIRINC_BOOST) -lboost_regex -o $@

.PHONY: clean
clean:
	rm -rf $(EXE)
