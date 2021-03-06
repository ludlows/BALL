# -*- Mode: C++; tab-width: 2; -*-
# vi: set ts=2:
#
# $Id: Makefile.in,v 1.40.4.2 2005/11/09 12:09:58 oliver Exp $
#

include config.mak

# 	library version number
LIB_VERSION	= 0

#	default library name
LIBNAME		= libBALL

# 	names of the static libraries
STATIC_LIBS	= libBALL.a 

#	names of the shared libraries
SHARED_LIBS	= $(STATIC_LIBS:.a=.so)

# names of the library object file lists
OBJECT_LISTS = $(STATIC_LIBS:.a=.objects)

.SUFFIXES:	.objects .a	.so

default:	lib

sharedlibs: $(SHARED_LIBS)

staticlibs: $(STATIC_LIBS)

install:	install.shared

BALLView:
	@( cd APPLICATIONS/BALLVIEW; make )

PB:
	@( cd APPLICATIONS/PB; make )

amber:
	@( cd APPLICATIONS/AMBER; make )

clean:
	@$(ECHO) "Cleaning up object files, cores,  and (uninstalled) library files..."
	@$(MAKE) subdirs SUBDIRS="$(DIRS) TEST BENCHMARKS APPLICATIONS" TARGET=default_clean
	@-$(RM) -f  $(STATIC_LIBS) $(SHARED_LIBS)
	@-$(RM) -rf .Dependencies */.Dependencies *.objects *~ 

lib_clean:
	@$(ECHO) "Removing all installed libraries..."
	@-$(RM) -rf ../lib


libBALL.objects:
	@$(ECHO) "Collecting object files for libraries..."
	@-rm -f *.objects
	@$(MAKE) subdirs SUBDIRS="$(DIRS)" TARGET=collectlib LIBNAME=$(LIBNAME)
	
libVIEW.objects:
	@$(ECHO) "Collecting object files for libraries..."
	@-rm -f *.objects
	@$(MAKE) subdirs SUBDIRS="$(DIRS)" TARGET=collectlib LIBNAME=$(LIBNAME)
	

collect:	$(OBJECT_LISTS)
	

libBALL.so: libBALL.objects
	$(DYNAR) $(ADD_DYNAROPTS_LIBBALL) $(DYNAROPTS) libBALL.so `cat libBALL.objects` $(PYTHON_LIBS) $(FFTW_LIBS) -lm

libVIEW.so: libVIEW.objects libBALL.so
	$(DYNAR) $(ADD_DYNAROPTS_LIBVIEW) $(DYNAROPTS) libVIEW.so `cat libVIEW.objects` -L. -lBALL $(QT_LIBOPTS) $(OPENGL_LIBOPTS) $(X11_LIBOPTS) $(PYTHON_LIBS) -lm

.objects.a:	
	@ $(ECHO) "Creating library $@ ..."
	@ $(AR) $(AROPTS) $*.a `cat $*.objects`
	@ echo $*.a

lib:	.Dependencies dirs collect $(SHARED_LIBS)
	@$(ECHO) " "
	@$(ECHO) "You may now install the libraries with the command"
	@$(ECHO) " "
	@$(ECHO) "   make install"
	@$(ECHO) " "
	@$(ECHO) "If you do not trust this Makefile, try"
	@$(ECHO) " "
	@$(ECHO) "   make -n install"
	@$(ECHO) " "
	@$(ECHO) "first to check where the libraries would be installed."
	@$(ECHO) " "
	@$(ECHO) "After you ran make install, please verify the validity"
	@$(ECHO) "of the installed library by executing"
	@$(ECHO) "   make test"
	@$(ECHO) " "
	@$(ECHO) "This will compile and run an extensive test suite."
	@$(ECHO) " "
	@$(ECHO) "If you use the shared libraries, don't forget to add"
	@$(ECHO) "the path "
	@$(ECHO) " $(BALL_PATH)/lib/$(BINFMT)"
	@$(ECHO) "to your LD_LIBRARY_PATH. "
	@$(ECHO) ""

install.static:	
	@$(ECHO) Installing static libraries under $(BALL_PATH)/lib/$(BINFMT)...
	@-mkdir $(BALL_PATH)/lib 2>/dev/null | cat >/dev/null
	@-mkdir $(BALL_PATH)/lib/$(BINFMT) 2>/dev/null | cat >/dev/null
	@$(MV) -f $(STATIC_LIBS) $(BALL_PATH)/lib/$(BINFMT)

install.shared:	
	@$(ECHO) Installing shared libraries under $(BALL_PATH)/lib/$(BINFMT)...
	@-mkdir $(BALL_PATH)/lib 2>/dev/null | cat >/dev/null 
	@-mkdir $(BALL_PATH)/lib/$(BINFMT) 2>/dev/null | cat >/dev/null
	@$(MV) -f $(SHARED_LIBS) $(BALL_PATH)/lib/$(BINFMT)

delete.static:	
	@$(ECHO) deleting static libraries under $(BALL_PATH)/lib/$(BINFMT)...
	@set -e; for i in $(STATIC_LIBS) ; do $(RM) $(BALL_PATH)/lib/$(BINFMT)/$$i 2>/dev/null ; done

delete.shared:	
	@$(ECHO) deleting shared libraries under $(BALL_PATH)/lib/$(BINFMT)...
	@set -e; for i in $(SHARED_LIBS) ; do $(RM) $(BALL_PATH)/lib/$(BINFMT)/$$i 2>/dev/null ; done

.Dependencies:
	@$(MAKE) depend

depend:
	@$(ECHO) "creating dependencies..."
	@$(ECHO) "" >.Dependencies
	@for i in $(DIRS) TEST BENCHMARKS APPLICATIONS ; do ( cd $$i; echo "" >.Dependencies; cd .. ) ; done
	@set -e; for i in $(DIRS) TEST BENCHMARKS APPLICATIONS ; do ( cd $$i; $(MAKE) depend; cd .. ) ; done


dirs:		dummy
	@$(MAKE) subdirs SUBDIRS="$(DIRS)" TARGET=default

subdirs:	dummy
	@set -e; for i in $(SUBDIRS); do cd $$i; if $(MAKE) $(TARGET); then cd .. ; else cd .. ; exit 1 ; fi ; done


test:	dummy
	@cd TEST ; $(MAKE) test

bench:	dummy
	@cd BENCHMARKS ; $(MAKE) bench


#
#    special targets - internal use only! 
#

# creates the documentation
doc:
	@config/makedoc

# creates the documentation for the website
doc-website:
	@config/makedoc -website

# creates a distributable .tgz file
dist:
	@config/makedistribution

# creates distribution for internal use - everything included
intdist:
	@config/makedistribution -int

# creates distribution for internal use - without documentation
intdistnodoc:
	@config/makedistribution -intnodoc

dummy:



configure:  config/configure.in config/configure.header
	@$(ECHO) "Creating new ./configure from config/configure.in"
	@autoconf config/configure.in > configure.tmp
	@cat configure.tmp|sed "2r config/configure.header" >configure
	@-rm configure.tmp
	@chmod 755 configure
 
