#
# BALL Python support Makefile
#
# $Id: Makefile,v 1.33.4.11 2005/11/12 14:16:29 oliver Exp $
#

DIRECTORY=PYTHON/EXTENSIONS

# allow processing of SIP files
.SUFFIXES: .sip

# default target: the BALL Python extensions in BALLmodule.so
all:	BALLmodule.so

#
include ../../common.mak
-include BALL.sbf

SIP_SOURCES=\
	$(sources)

SIP_OBJECTS	=	\
	$(sources:.C=.o)

SIP_HEADERS = \
	$(PYBALL_SOURCES:.C=.h) \

SIP_HEADER_FILES =\
	$(headers)
	
SIP_DIR=sip_files
SIP_OPTS=-c $(SIP_DIR) -e -s ".C" -b BALL.sbf
SIPFLAGS += $(PYTHON_INCLUDES) $(SIP_INCLUDES)
CXXFLAGS += $(SIPFLAGS)

BALL_INCLUDES+=$(VIEW_INCLUDES)
BALL_LIBS=$(VIEW_LIBS)
ADD_INCLUDES=$(PYTHON_INCLUDES) $(SIP_INCLUDES)
 
PYBALL_SOURCES = \
	extractors.C \
	pyBALLSipHelper.C\
	pyCompositeList.C\
	pyRepresentationList.C\
	pyIndexList.C\
	pyAtomList.C \
	pyBondList.C \
	pyConstant.C \
	pyMoleculeList.C \
	pyNucleotideList.C \
	pyProteinList.C \
	pySecondaryStructureList.C \
	pyAtomContainerList.C \
	pyChainList.C \
	pyFragmentList.C \
	pyNucleicAcidList.C \
	pyPDBAtomList.C \
	pyResidueList.C \
	pyGeometricObjectList.C
	
#pyColorRGBAVector.C

PYBALL_OBJECTS  = $(PYBALL_SOURCES:.C=.o)
 
CPP_SOURCES = $(SIP_SOURCES:.h=.C) $(PYBALL_SOURCES)
CPP_OBJECTS = $(SIP_OBJECTS) $(PYBALL_OBJECTS)

BALLmodule.so: sipBALLcmodule.o

sip:
	@if test ! -d $(SIP_DIR) ; then mkdir $(SIP_DIR) 2>/dev/null ; fi
	@echo "running SIP..."
	@$(SIP) $(SIP_OPTS) BALL.sip
	@echo "removing SIP headers..."
	@-rm -f $(SIP_DIR)/tmp 2>/dev/null
	@for i in $(SIP_DIR)/*; do mv $$i $(SIP_DIR)/tmp ; cat $(SIP_DIR)/tmp |sed 4d > $$i ; rm -rf $(SIP_DIR)/tmp 2>/dev/null ; done
	@echo "updating modified files..."
	@for i in $(SIP_DIR)/*; do if cmp -s $$i `basename $$i`; then true ; else echo " - $$i"; cp $$i . ; fi ; done


BALLmodule.so:  $(CPP_OBJECTS)
	$(DYNAR) $(DYNAROPTS_PYMODULE) BALLmodule.so $(SIP_OBJECTS) $(PYBALL_OBJECTS) $(SIP_LIB) -L$(BALL_PATH)/lib/$(BINFMT) $(BALL_LIBS) $(PYTHON_LIBS)

install:	BALLmodule.so
	cp BALLmodule.so  $(BALL_PATH)/lib/$(BINFMT)
	 
CLEAN_FILES = \
	BALLmodule.so\
	$(CPP_OBJECTS) \
	$(SIP_SOURCES) \
	$(SIP_HEADERS) \
	$(SIP_HEADER_FILES) \
	$(SIP_OBJECTS)
 
clean:
	@echo "removing object files, libraries, and source/header files produced by SIP..."
	@-$(RM) -f $(CLEAN_FILES) 2>/dev/null
	@-$(RM) -fr $(SIP_DIR)

depend: sip  Makefile.inc .Dependencies

-include .Dependencies
