# 
# $Id: common.mak.in,v 1.14 2004/02/17 14:46:04 oliver Exp $
#

include /home/wmcs/work/BALL/source/config.mak

SOURCES		= $(CPP_SOURCES) $(MOC_SOURCES) $(ADD_CPP_SOURCES) $(ADD_MOC_SOURCES)
OBJECTS 	= $(SOURCES:.C=.o)

# support for conditional Pythone source compilation
#  PYTHON_SUPPORT has to be defined as empty (do not build Python stuff) 
#  or as true (build everything)
# Makefiles have to define PYTHON_MOC_SOURCES and PYTHON_CPP_SOURCES for 
# for this purpose
ADD_CPP_SOURCES = $(PYTHON_SUPPORT:true=$(PYTHON_CPP_SOURCES))
ADD_MOC_SOURCES = $(PYTHON_SUPPORT:true=$(PYTHON_MOC_SOURCES))

default:	echodir $(OBJECTS)

echodir:
	@echo "entering $(DIRECTORY)..."


# .dh -- dummy header files for 
# .il -- Intel CC otimizer files
default_clean:
	@-$(RM) -rf *.dh *il *~  core core.* $(OBJECTS) $(MOC_SOURCES) $(UIC_SOURCES) $(ADD_MOC_SOURCES) $(EXECUTABLES) $(EXECUTABLE)


collectlib:	default
	@echo "collecting objects from $(DIRECTORY)..."
	@for i in $(OBJECTS); do echo $(DIRECTORY)/$$i >>/home/wmcs/work/BALL/source/$(THISLIB).objects ; done


depend: $(CPP_SOURCES)
	@echo "creating dependencies for $(DIRECTORY)..."
	@echo "" > .Dependencies
	@/usr/bin/g++ -M  -pipe -fPIC   $(CPP_SOURCES) $(ADD_CPP_SOURCES) $(STD_CPP_INCLUDES) $(ADD_INCLUDES) $(BALL_INCLUDES)  >.Dependencies
