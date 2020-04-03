CC = g++ 
PROFILE =
LAPACK =
CHECKMKL = $(shell gcc -lmkl_rt 2>&1 | grep "cannot find")
CHECKLAPACK = $(shell gcc -llapack 2>&1 | grep "cannot find")
ifeq ($(LAPACK),)
  ifeq ($(CHECKMKL),)
    LAPACK = -lmkl_rt
  endif
endif
ifeq ($(LAPACK),)
  ifeq ($(CHECKLAPACK),)
    LAPACK = -lblas -llapack
  endif
endif
#PROFILE = -pg
#PROFILE = -g
CFLAGS := -c -Wall -Wextra -pedantic -std=gnu++11 -Ofast $(PROFILE)
LDFLAGS = $(PROFILE)
#comment out to deactivate debug and asserts
#CFLAGS := $(CFLAGS) -D NDEBUG
#use rational numbers from boost
#CFLAGS := $(CFLAGS) -D _RATIONAL
#use lapack
ifdef LAPACK
  CFLAGS := $(CFLAGS) -D _LAPACK
endif
INCLUDES=$(LAPACK)
# program name
MAIN = dumpham
# OS type
UNAME_S := $(shell uname -s)
# base directory
BASE=$(PWD)
ifeq ($(UNAME_S),Linux)
  BASE=$(shell dirname $$(readlink -f Makefile))
endif
ifeq ($(findstring CYGWIN, $(UNAME_S)),CYGWIN)
  BASE=$(shell dirname $$(readlink -f Makefile))
endif
ifeq ($(UNAME_S),Darwin)
  # uses greadlink from coreutils
  BASE=$(shell dirname $$(greadlink -f Makefile))
endif
in=input
out=FCIDUMP
# files to be linked to working-directory
FILIN=fcidumpcalc
DIR = src
OBJ0 = main.o FCIdump.o hdump.o refdet.o odump.o integs.o finput.o inpline.o utilities.o globals.o density_mat.o fock_mat.o
OBJ = $(patsubst %,$(DIR)/%,$(OBJ0))
SRC = $(OBJ:.o=.cpp)

all : $(MAIN)

$(MAIN) : $(OBJ)
	$(CC) $(LDFLAGS) $(OBJ) $(INCLUDES) -o $(MAIN)

%.o : %.cpp 
	$(CC) $(CFLAGS) $< -o $@

clean : 
	rm -rf $(MAIN) $(OBJ) 
veryclean :
	git clean -dfx  
mol :
	 $(MAIN) $(in).dh
	 fcidumpcalc $(out)

base : $(FILIN)
ifneq ($(BASE),$(PWD))
	 @test -e $(MAIN) || ln -s $(BASE)/$(MAIN) .
endif

$(FILIN) :
ifneq ($(BASE),$(PWD))
	 @test -e $@ || ln -s $(BASE)/$@ .
endif

depend: 
	 makedepend -- $(INCLUDES) $(CFLAGS) -Y -- $(SRC) 

# DO NOT DELETE THIS LINE -- make depend needs it

src/main.o: src/argpars.h src/utilities.h src/globals.h src/hdtypes.h
src/main.o: src/finput.h src/inpline.h src/hdump.h src/FCIdump.h src/pgsym.h
src/main.o: src/integs.h src/refdet.h src/odump.h src/density_mat.h
src/main.o: src/fock_mat.h
src/FCIdump.o: src/FCIdump.h
src/hdump.o: src/hdump.h src/globals.h src/hdtypes.h src/utilities.h
src/hdump.o: src/FCIdump.h src/pgsym.h src/integs.h src/refdet.h
src/refdet.o: src/refdet.h src/globals.h src/hdtypes.h src/utilities.h
src/refdet.o: src/pgsym.h
src/odump.o: src/odump.h src/pgsym.h src/globals.h src/hdtypes.h
src/odump.o: src/utilities.h src/integs.h src/refdet.h src/inpline.h
src/integs.o: src/integs.h src/pgsym.h src/globals.h src/hdtypes.h
src/integs.o: src/utilities.h
src/finput.o: src/finput.h src/utilities.h src/globals.h src/hdtypes.h
src/finput.o: src/inpline.h src/hdump.h src/FCIdump.h src/pgsym.h
src/finput.o: src/integs.h src/refdet.h src/odump.h src/density_mat.h
src/finput.o: src/fock_mat.h
src/inpline.o: src/inpline.h src/utilities.h src/globals.h src/hdtypes.h
src/utilities.o: src/utilities.h src/globals.h src/hdtypes.h
src/globals.o: src/globals.h src/hdtypes.h
src/density_mat.o: src/density_mat.h src/globals.h src/hdtypes.h
src/density_mat.o: src/utilities.h src/odump.h src/pgsym.h src/integs.h
src/density_mat.o: src/refdet.h src/inpline.h
src/fock_mat.o: src/fock_mat.h src/integs.h src/pgsym.h src/globals.h
src/fock_mat.o: src/hdtypes.h src/utilities.h src/density_mat.h src/odump.h
src/fock_mat.o: src/refdet.h src/inpline.h src/hdump.h src/FCIdump.h
