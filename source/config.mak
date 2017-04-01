# This is config.mak.in/config.mak
# Please do ONLY change config.mak.in, as config.mak is automagically 
# created by configure from config.mak.in
#
# $Id: config.mak.in,v 1.35.2.3 2005/11/09 12:09:58 oliver Exp $
#

# First, some basic unix commands, including the C++ compiler
#

# archiver ar (to create static libraries)
AR	=/usr/bin/ar
RANLIB = echo

# dynamic archiver (to create shared libraries)
DYNAR = /usr/bin/g++

# options for library generation
AROPTS  = cru

# options for shared library generation
DYNAROPTS =  -shared -fPIC -o
DYNAROPTS_PYMODULE =  -shared -fPIC -o
ADD_DYNAROPTS_LIBBALL = 
ADD_DYNAROPTS_LIBVIEW = 

# suffix of shared library filenames
SHARED_LIB_SUFFIX = so
EXEEXT = 
 
# C++ compiler
CXX=/usr/bin/g++
# unix command rm, used to remove obsolete files (target clean)
RM=/bin/rm
# unix command cp, used to copy files
CP=/bin/cp
# unix command mv, used to install libraries (target install.*)
MV=/bin/mv
# unix command strip, used to strip away unused stuff from object files
STRIP=@STRIP@
# the echo command (on some platforms with an additional -e to expand '\n')
ECHO=echo
# flex and bison
LEX	=flex
YACC=bison -y

OS	     =Linux
# name of your OS
OSREV	   =4.4.0-64-generic
# full revision number
OSMAJOR	 =4
# major revision number (full revision number cut at the first dot)
BINFMT   =Linux-x86_64-g++_4.8.4
# binary format type (platform-architecture-compiler)

# BALL_PATH is the BALL root directory (the directory in which include and source reside)
# BALL_INCLUDES contains compiler flags to set the BALL include paths
# and BALL_LIBS is set to $BALL_PATH/lib/$BINFMT
#
BALL_PATH		    = /home/wmcs/work/BALL
BALL_INCLUDES		= -I/home/wmcs/work/BALL/include 
BALL_LIBS       =  -L/home/wmcs/work/BALL/lib/Linux-x86_64-g++_4.8.4 -lBALL  -lm

# USE_VIEW contains true if VIEW support is to be built, false otherwise
# VIEW_PLATFORM either contains Mesa or OpenGL (so far)
#
USE_VIEW		    = false
VIEW_PLATFORM 	= OpenGL
VIEW_INCLUDES	  = 
VIEW_LIBS		    =      -lm
BALL_HAS_XDR				= true
BALL_HAS_VIEW				= 
BALL_HAS_FFTW       = 
BALL_HAS_FLEX_YACC	= true
OPENGL_LIBPATH      = 
OPENGL_LIBOPTS      = 
QT_LIBPATH          = 
QT_MT_SUFFIX        = -mt
QT_LIBOPTS          = 
MOC                 = moc
UIC                 = uic
X11_LIBPATH         = 
X11_LIBOPTS         = 
X11_LIBS            = 
SIP            			= 
SIP_LIB        			= 
SIP_INCLUDES   			= 
PYTHON_INCLUDES			= 
PYTHON_SUPPORT      =  
PYTHON_LIBS					= 
PYTHON_EXECUTABLE	  =	
FFTW_LIBS						= 

#
# Support for valgrind, a memory leak checker.
# Try "make valgrind" in source/TEST.
#
VALGRIND=valgrind
VALGRIND_OPTS=-v --leak-check=yes --leak-resolution=high

# These are standard flags. CXXFLAGS to compile, LDFLAGS to link programs
# DEFINES and DEBUG are both part of CXXFLAGS
# ADD_CXXFLAGS are also passed on to the C++ compiler
# but are NOT passed on to makedepend
# This is the place for compiler options that are not understood
# by makedepend (e.g. -Olimit using IRIX DCC)
#
CXXFLAGS		    =  -pipe -fPIC  
LDFLAGS			    = 
ADD_CXXFLAGS		= 

# all the libraries in linker usable format (-l..)
LIBS			=  -L/home/wmcs/work/BALL/lib/Linux-x86_64-g++_4.8.4 -lBALL  -lm
#
DEFINES   = 
#
CXXFLAGS_D		=  -Wall -W -pedantic -Wno-long-long
CXXFLAGS_O		=  -O3 -Wall -W -pedantic -Wno-long-long

#
# additional CXXFLAGS for building library objects
LIB_CXXFLAGS = 

#
# additional CXXFLAGS for building applications, tests, examples
# used instead of LIB_CXXFLAGS in APPLICATIONS, EXAMPLES, TEST, BENCHMARKS
NONLIB_CXXFLAGS = 

# MODE_FLAGS is set to DEBUG_FLAGS or OPTIMIZE_FLAGS, this
# can be controlled by specifying the configure switch 
# --enable-optimization
CPP_MODE_FLAGS	=  -O3 -Wall -W -pedantic -Wno-long-long

# CPP_MODE_FLAGS_NO_OPTIMIZATION is used in the (luckily quite rare!)
# cases where the compiler crashes when performing optimization.
# In this case, the source files will be compiled using CPP_MODE_FLAGS_NO_OPTIMIZATION
# instead of CPP_MODE_FLAGS
CPP_MODE_FLAGS_NO_OPTIMIZATION = 

# if g++ is used, configure tries to determine its standard include directories
# this is used to prevent trouble when calling makedepend. Otherwise makedepend
# will produce loads of (irrelevant) warnings
STD_CPP_INCLUDES	= 

#  These are the major source directories.
#  They are entered in this order and "make all" is executed.
#  To build the library each of these directories is entered und "make libadd" is executed
#  to add the subdirectories' objects to the library.
#
DIRS	= COMMON CONCEPT DATATYPE ENERGY FORMAT KERNEL MATHS MOLMEC NMR QSAR SOLVATION STRUCTURE SYSTEM  

#
.SUFFIXES: .C .y .l .o .vgr _moc.C .ui .h Data.C Data_moc.C Data.h .dh Data.C

#
# Rule to create .o files from .C files
#
.C.o:
	$(CXX) $(CXXFLAGS) $(ADD_CXXFLAGS) $(CPP_MODE_FLAGS) $(LIB_CXXFLAGS) $(BALL_INCLUDES) $(ADD_INCLUDES) -c $*.C -o $@

#
# Rule to create .C files from the BISON .y files
#
.y.C:
	$(YACC) -p $(PARSER_PREFIX) -d $*.y && $(MV) y.tab.h $*.h && $(MV) y.tab.c $@

#
# Rule to create .C files from the FLEX .l files
#
.l.C:
	$(LEX) -P$(PARSER_PREFIX) -o$@ $*.l

#
# Rule to create header files from .ui files (using UIC)
# .dh files are dummy files required due to the stupidity 
# of some implementations of make (e.g. SUN)
# They are just empty files that exist whenever the corresponding
# *Data.h files exists.
#
.ui.dh:
	$(UIC) -o $(UIC_DIR)/$*Data.h $< && touch $*.dh

.dhData.C:
	$(UIC) -o $@ -impl $(UIC_DIR)/$*Data.h $*.ui

#
# Rule to create moc files from header files (using MOC)
#
.dhData_moc.C: $*.dh
	$(MOC) -o $@ $(UIC_DIR)/$*Data.h

.dh_moc.C:  $*.dh
	$(MOC) -o $@ $(HEADER_DIR)/$*.h

#
# Rule to create valgrind output from compiled test
#
.vgr:	$*
	$(VALGRIND) $(VALGRIND_OPTS) $*
