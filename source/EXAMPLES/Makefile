# -*- Mode: C++; tab-width: 2; -*-
# vi: set ts=2:
#
# $Id: Makefile,v 1.21.2.2 2005/11/01 08:52:51 oliver Exp $
#

#========================
# BALL/EXAMPLES Makefile
#========================


include ../common.mak
-include .Dependencies

DIRECTORY=EXAMPLES
EXECUTABLES=proteinMapper PDB2discover amberEnergy bindingEnergy ramachandran PDBChecker addHydrogens SASAreas energy
VIEW_EXECUTABLES=$(BALL_HAS_VIEW:true=BALLVIEWClient)

BALL_INCLUDES=$(VIEW_INCLUDES) $(PYTHON_INCLUDES)


default:	$(EXECUTABLES) $(VIEW_EXECUTABLES)

# remove the library specific flags
LIB_CXXFLAGS = $(NONLIB_CXXFLAGS)

clean:	default_clean
	@-$(RM) $(EXECUTABLES) *.o 2>/dev/null

proteinMapper:	proteinMapper.o 
	$(CXX)	$(LDFLAGS) -o $@ $@.o $(LIBS)

PDB2discover:	PDB2discover.o 
	$(CXX)	$(LDFLAGS) -o $@ $@.o $(LIBS)

amberEnergy:	amberEnergy.o 
	$(CXX)	$(LDFLAGS) -o $@ $@.o $(LIBS)

energy:	energy.o 
	$(CXX)	$(LDFLAGS) -o $@ $@.o $(LIBS)

bindingEnergy:	bindingEnergy.o 
	$(CXX)	$(LDFLAGS) -o $@ $@.o $(LIBS)

BALLVIEWClient:	BALLVIEWClient.o 
	$(CXX)	$(LDFLAGS) -o $@ $@.o $(VIEW_LIBS)

ramachandran:	ramachandran.o 
	$(CXX)	$(LDFLAGS) -o $@ $@.o $(LIBS)

PDBChecker:	PDBChecker.o 
	$(CXX)	$(LDFLAGS) -o $@ $@.o $(LIBS)

addHydrogens:	addHydrogens.o
	$(CXX)	$(LDFLAGS) -o $@ $@.o $(LIBS)

SASAreas:	SASAreas.o 
	$(CXX)	$(LDFLAGS) -o $@ $@.o $(LIBS)
