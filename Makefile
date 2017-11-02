.EXPORT_ALL_VARIABLES:

.PHONY: clean all

LIB_DIR = $(HOME)/lib

ROOTCFLAGS   := $(shell root-config --cflags)
ROOTLIBS     := $(shell root-config --libs)
ROOTGLIBS    := $(shell root-config --glibs)
ROOTINC      := -I$(shell root-config --incdir)

COMMON_DIR = $(HOME)/TINAanalysis/common


BASELIBS  = -lm $(ROOTLIBS) $(ROOTGLIBS) -L$(LIB_DIR)
ALLIBS  = $(BASELIBS) -lCommandLineInterface -lKinematics

CPP             = g++
#CFLAGS	        = -g -O3 $(ROOTCFLAGS)
CFLAGS		= -pedantic -Wall -Wno-long-long -g -O3 $(ROOTCFLAGS) -fPIC
DFLAGS		= -Wall -Wno-long-long -g -O3 $(ROOTCFLAGS) -fPIC

INCLUDES        = -I./ -I$(COMMON_DIR) 
LFLAGS		= -g -fPIC
LIBS 		= $(ALLIBS)

CFLAGS += -Wl,--no-as-needed
LFLAGS += -Wl,--no-as-needed
DFLAGS += -Wl,--no-as-needed

O_FILES = Nucleus.o \
	Compound.o \
	Reconstruction.o 

all: ELoss
	echo Done

ELoss: ELoss.cc $(O_FILES)
	$(CPP) $(CFLAGS) $(INCLUDES)  $^ $(LIBS) -o $@
	cp ELoss $(HOME)/bin

lib%.so: Nucleus.o NucleusDictionary.o Compound.o CompoundDictionary.o Reconstruction.o ReconstructionDictionary.o
	$(CPP) $(LFLAGS) -shared -Wl,-soname,$@ -o $(LIB_DIR)/$@ $^ $(BASELIBS) -lc

%.o: %.cc %.hh
	@echo Default .o rule
	$(CPP) $(CFLAGS) $(INCLUDES) -c $< $(LIBS) -o $@

NucleusDictionary.o: NucleusDictionary.cc NucleusDictionary.h
	 $(CPP) -p -fPIC $(DFLAGS) -c $<

NucleusDictionary.cc: Nucleus.hh NucleusLinkDef.h 
	 rm -f NucleusDictionary.cc NucleusDictionary.h; rootcint -f $@ -c Nucleus.hh NucleusLinkDef.h 

CompoundDictionary.o: CompoundDictionary.cc CompoundDictionary.h
	 $(CPP) -p -fPIC $(DFLAGS) -c $<

CompoundDictionary.cc: Compound.hh CompoundLinkDef.h 
	 rm -f CompoundDictionary.cc CompoundDictionary.h; rootcint -f $@ -c Compound.hh CompoundLinkDef.h 

ReconstructionDictionary.o: ReconstructionDictionary.cc ReconstructionDictionary.h
	 $(CPP) -p -fPIC $(DFLAGS) -c $<

ReconstructionDictionary.cc: Reconstruction.hh ReconstructionLinkDef.h 
	 rm -f ReconstructionDictionary.cc ReconstructionDictionary.h; rootcint -f $@ -c Reconstruction.hh ReconstructionLinkDef.h 

clean:
	rm *.o *Dictionary.* ELoss
