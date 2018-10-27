CC = g++ 
PROFILE =
#PROFILE = -pg
#PROFILE = -g
CFLAGS := -c -Wall -Wextra -pedantic -std=gnu++11 -Ofast $(PROFILE)
LDFLAGS = $(PROFILE)
#comment out to deactivate debug and asserts
#CFLAGS := $(CFLAGS) -D NDEBUG
#use rational numbers from boost
#CFLAGS := $(CFLAGS) -D _RATIONAL
INCLUDES=
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
OBJ0 = main.o FCIdump.o hdump.o odump.o finput.o inpline.o utilities.o globals.o
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

src/main.o: src/argpars.h src/utilities.h src/globals.h src/finput.h
src/main.o: src/inpline.h src/hdump.h src/FCIdump.h src/odump.h
src/FCIdump.o: src/FCIdump.h
src/hdump.o: src/hdump.h src/globals.h src/utilities.h src/FCIdump.h
src/odump.o: src/odump.h src/globals.h src/utilities.h src/inpline.h
src/finput.o: src/finput.h src/utilities.h src/globals.h src/inpline.h
src/finput.o: src/hdump.h src/FCIdump.h src/odump.h
src/inpline.o: src/inpline.h src/utilities.h src/globals.h
src/utilities.o: src/utilities.h src/globals.h
src/globals.o: src/globals.h
