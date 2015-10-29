.EXPORT_ALL_VARIABLES:

.PHONY: clean all

LIB_DIR = $(HOME)/lib

ROOTCFLAGS   := $(shell root-config --cflags)
ROOTLIBS     := $(shell root-config --libs)
ROOTGLIBS    := $(shell root-config --glibs)
ROOTINC      := -I$(shell root-config --incdir)

COMMON_DIR = $(HOME)/common


ALLIBS  = -lm $(ROOTLIBS) $(ROOTGLIBS) -L$(LIB_DIR) -lCommandLineInterface 
CC		= gcc
CPP             = g++
CFLAGS		= -g -O3 $(ROOTCFLAGS)

INCLUDES        = -I./ -I$(COMMON_DIR)
LFLAGS		= -g 
LIBS 		= $(ALLIBS)

O_FILES = Reconstruction.o \
	Nucleus.o \
	Compound.o 

LIBRARIES = $(LIB_DIR)/libCommandLineInterface.so 

all: ELoss
	echo Done

ELoss: ELoss.cc $(O_FILES)
	$(CPP) $(CFLAGS) $(INCLUDES)  $^ $(LIBS) -o $@
	cp ELoss $(HOME)/bin

%.o: %.cc %.hh
	@echo Default .o rule
	$(CPP) $(CFLAGS) $(INCLUDES) -c $< -o $@

