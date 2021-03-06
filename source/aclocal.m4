dnl -*- Mode: C++; tab-width: 1; -*-
dnl vi: set ts=2:
dnl
dnl		$Id: aclocal.m4,v 1.69.2.24 2006/01/09 22:17:03 oliver Exp $
dnl
dnl		Autoconf M4 macros used by configure.ac.
dnl

dnl
dnl Determine the current version of the code and print it
dnl to the console to simplify debugging of submitted configure
dnl output.
dnl
AC_DEFUN(CF_VERSION_CHECK,[
MACHINE=`uname -a`
AC_MSG_RESULT(This is PROJECT $PROJECT[]_VERSION_STRING compiling.)
AC_MSG_RESULT(Host: $MACHINE)
])

dnl
dnl		display the license and abort if not accepted
dnl		we create the file config.lic if the license was
dnl		accepted and do not show the license the second time
dnl

AC_DEFUN(CF_CHECK_LICENSE,[
AC_PATH_PROG(PAGER, more, no)
if test "${PAGER}" = "no" ; then
	PAGER=cat
fi
if test ! -f config.lic ; then
	${PAGER} COPYRIGHT
	echo " "
	echo "Do you accept these terms (y/n)?"
	answer=""
	while test "$answer" != "y" -a "$answer" != "n" ; do
		read answer 
	done
	if test "$answer" = "n" ; then	
		exit
	else
		echo "accepted" > config.lic
	fi
fi
])

dnl		define a macro to remove the directory name
dnl		from a fully specified path
dnl
AC_DEFUN(CF_BASENAME,[
TMP__DIR="$1/"
while test "${TMP__DIR}" != "" ; do
	TMP__NAME=`echo ${TMP__DIR}|${CUT} -d/ -f1`
	TMP__DIR=`echo ${TMP__DIR}|${CT} -d/ -f2-`
	if test "${TMP__DIR}" = "${TMP__NAME}" ; then
		TMP__DIR=""
	fi
	done
])

dnl    define a macro to abort configure, print an appropriate error message
dnl    and package up the current stuff relevant to diagnosis into a tar
dnl    file.
AC_DEFUN(CF_ERROR,[
AC_MSG_RESULT()
AC_MSG_RESULT([Configure failed. If you cannot solve your problem with the aid])
AC_MSG_RESULT([of the above error message, please contact the ]PROJECT[ mailing list])
AC_MSG_RESULT([or the ]PROJECT[ developers. Please enclose the file 'conf.diag.tar'])
AC_MSG_RESULT([which has been created in source. It contains the relevant])
AC_MSG_RESULT([files from this configure run. In most cases, the information])
AC_MSG_RESULT([is necessary to diagnose what went wrong. This file contains])
AC_MSG_RESULT([information about your system setup and versions of compilers])
AC_MSG_RESULT([and other tools installed in your system.])
AC_MSG_RESULT()
CF_CONF_DIAG_TAR
AC_ERROR(Aborted.)
])

AC_DEFUN(CF_CONF_DIAG_TAR,[
TARFILE=conf.diag.tar
if test -f $TARFILE ; then 
	${RM} $TARFILE ; 
fi
FILES="configure configure.ac aclocal.m4 config.log"
if test -f config.h ; then FILES="$FILES config.h" ; fi
if test -f config.mak ; then FILES="$FILES config.mak" ; fi
if test -f common.mak ; then FILES="$FILES common.mak" ; fi
tar cf $TARFILE $FILES
])

dnl    define a macro to inform the user about failed tests for programs
dnl    it checks for the unix command given as second parameter and
dnl    sets the shell variable given as second parameter to its absolute path
dnl 
AC_DEFUN(CF_MSG_PATH_PROG,[
AC_PATH_PROG($1,$2,no)
if test $$1 = no ; then
	AC_MSG_RESULT()
	AC_MSG_RESULT([This script requires the unix command $2, but cannot find it.])
	AC_MSG_RESULT([Please add the correct path to $2 to your \$PATH variable])
	AC_MSG_RESULT([and restart configure.])
	CF_ERROR
fi
])

dnl
dnl    define macro to search for header files that may be somewhere in the filesystem
dnl    if ${FIND}!=- (i.e. it has been set BEFORE py AC_PATH_PROG!) find will be used
dnl    too if everything fails - this may take some time...
dnl
dnl    syntax: AC_FIND_HEADER(<PATH_VAR>,<header.h>,<additional dirnames>)
dnl    
dnl        PATH_VAR will be set to the include path or to empty string (if not found)
dnl        header.h is the header file name (e.g. wait.h, GL/gl.h)
dnl        additional dirnames are included in searches these should be absolute names.
dnl 

AC_DEFUN(CF_FIND_HEADER,[
_INCLUDES=

dnl    immediate return on predefined directory (read from file?)
if test "${$1}" != "" ; then
_INCLUDES=${$1}
fi

if test "${_INCLUDES}" = "" ; then
for i in /usr/include /opt/include $3 ; do
if test -f "$i/$2" -a "${_INCLUDES}" = ""; then
_INCLUDES="$i"
fi
done
fi

if test "${_INCLUDES}" = "" ; then
for i in /usr/*/include /opt/*/include ; do
if test -f "$i/$2" -a "${_INCLUDES}" = ""; then
_INCLUDES="$i"
fi
done
fi

if test "${_INCLUDES}" = "" ; then
for i in /opt/*/*/include /usr/*/*/include /usr/local/*/*/include ; do
if test -f "$i/$2" -a "${_INCLUDES}" = ""; then
_INCLUDES="$i"
fi
done
fi

if test "${_INCLUDES}" = "" -a "${FIND}" != "-" ; then
if test "${FIND_KNOWS_PATH}" = false ; then
FIND_OPT="-name"
_TMP_FIND_NAME="$2"
while test _`${EGREP} / $_TMP_FIND_NAME`_ != __ ; do
_TMP_FIND_NAME=`echo ${_TMP_FIND_NAME}|${CUT} -d/ -f2-`
done

FIND_ARG="\*${_TMP_FIND_NAME}\*"
else
FIND_OPT="-path"
FIND_ARG="\*$2\*"
fi

_TMP=`${FIND} /usr ${FIND_OPT} ${FIND_ARG} -print 2>/dev/null`
for j in ${_TMP} ; do
if test "${_INCLUDES}" = "" ; then
_INCLUDES=`echo $j|${SED} "s/include\/.\*\$/include/"`
fi
done

if test "${_INCLUDES}" = "" ; then
_TMP=`${FIND} /opt ${FIND_OPT} ${FIND_ARG} -print 2>/dev/null`
for j in ${_TMP} ; do
if test "${_INCLUDES}" = "" ; then
	_INCLUDES=`echo $j|${SED} "s/include\/.\*\$/include/"`
fi
done
fi

if test "${_INCLUDES}" = "" -a "$3" != ""; then
for i in $3 /dev/null ; do
_TMP=`${FIND} $i ${FIND_OPT} ${FIND_ARG} -print 2>/dev/null`
for j in ${_TMP} ; do
	if test "${_INCLUDES}" = "" ; then
		_INCLUDES=`echo $j|${SED} "s/include\/.\*\$/include/"`
	fi
done
done
fi
fi

$1="${_INCLUDES}"
])


dnl
dnl    define macro to search for libraries that may be somewhere in the filesystem
dnl    if ${FIND}!=- (i.e. it has been set BEFORE py AC_PATH_PROG!) find will be used
dnl    too if everything fails - this may take some time...
dnl
dnl    syntax: CF_FIND_LIB(<PATH_VAR>,<libXXX>,<additional dirnames>)
dnl    
dnl        PATH_VAR will be set to the library path or to empty string (if not found)
dnl        libXXX is the header file name (e.g. libGLUT, libGL) .a, .so etc. should be omitted
dnl        additional dirnames are included in searches these should be absolute names.
dnl 

AC_DEFUN(CF_FIND_LIB,[
_LIBS=

dnl   immediate "return" on preset directory (read from file?)
if test "${$1}" != "" ; then
_LIBS=${$1}
fi

if test "$3" = "" ; then
_LIBDIRS="/usr/lib /opt/lib /usr/local/lib"
else
_LIBDIRS="$3"
fi

if test "${_LIBS}" = "" ; then
for i in ${_LIBDIRS} ; do
for j in $i/$2.* ; do
if test -f "$j" -a "${_LIBS}" = ""; then
	_LIBS="$i"
fi
done
done
fi

if test "${_LIBS}" = "" ; then
for i in /usr/lib /opt/lib /usr/*/lib /opt/*/lib ; do
for j in $i/$2.* ; do
if test -f "$j" -a "${_LIBS}" = ""; then
	_LIBS="$i"
fi
done
done
fi

if test "${_LIBS}" = "" ; then
for i in /opt/*/*/lib /usr/*/*/lib /usr/local/*/*/lib ; do
for j in $i/$2.* ; do
if test -f "$j" -a "${_LIBS}" = ""; then
	_LIBS="$i"
fi
done
done
fi

if test "${_LIBS}" = "" -a "${FIND}" != "-" ; then
if test "${_LIBS}" = "" -a "$3" != ""; then
for i in $3 /dev/null; do
if test "${_LIBS}" = "" ; then
	_TMP=`${FIND} $i -name "$2.*" -print 2>/dev/null`
	for j in ${_TMP} ; do
		if test "${_LIBS}" = "" ; then
			_LIBS=`echo $j|${SED} "s/\/$2\\.*/\//"`
		fi
	done
fi
done
fi

if test "${_LIBS}" = "" ; then
_TMP=`${FIND} /opt -name "$2.*" -print 2>/dev/null`
for j in ${_TMP} ; do
if test "${_LIBS}" = "" ; then
	_LIBS=`echo $j|${SED} "s/\/$2\\.*/\//"`
fi
done
fi

if test "${_LIBS}" = "" ; then
_TMP=`${FIND} /usr -name "$2.*" -print 2>/dev/null`
for j in ${_TMP} ; do
if test "${_LIBS}" = "" ; then
	_LIBS=`echo $j|${SED} "s/\/$2\\.*/\//"`
fi
done
fi

fi

$1="${_LIBS}"
])



dnl
dnl   check whether "echo" understands "-n" (required on some
dnl   platforms to expand '\n')
dnl
AC_DEFUN(CF_CHECK_ECHO,[
AC_MSG_CHECKING(whether echo accepts -e)
if `/bin/sh -c "echo -e \"\n\"" >/dev/null 2>&1` ; then
if test "`/bin/sh -c echo -e 2>&1`" = "" -a "`/bin/sh -c echo -e OK`" = "OK" ; then
ECHO_COMMAND="echo -e"
AC_MSG_RESULT(yes)	
else
ECHO_COMMAND="echo"
AC_MSG_RESULT(no)	
fi
else
ECHO_COMMAND="echo"
AC_MSG_RESULT(no)
fi
AC_SUBST(ECHO_COMMAND)
])


dnl
dnl   check whether find can be called with the parameter -path
dnl   (needed to find headers in a certain path like GL/libgl.h
dnl
AC_DEFUN(CF_CHECK_FIND,[
if test "${FIND}" != "no" ; then
RESULT=`${FIND} KERNEL -path . -print 2>&1`
if test "${RESULT}" != "" ; then     dnl    did get an error message ... bad.
FIND_KNOWS_PATH=false
else
FIND_KNOWS_PATH=true
fi
fi
])

dnl
dnl    determine OS and architecture and all this stuff
dnl
AC_DEFUN(CF_DETECT_OS,[
AC_SUBST(OSMAJOR)
AC_SUBST(OS)
AC_SUBST(OSREV)
AC_SUBST(BINFMT)
AC_SUBST(BINFMT_PATH)
AC_SUBST(ARCHITECTURE)
RANLIB="echo"

AC_MSG_CHECKING(your OS)
OS=`${UNAME} -s`
OSREV=`${UNAME} -r`
OSMAJOR=`echo $OSREV|${CUT} -d"." -f1`

dnl		default...
BINFMT="${OS}"

if test "$OS" = SunOS ; then
if test "$OSMAJOR" = 5 ; then
OS=Solaris
ARCHITECTURE=`${UNAME} -p`
BINFMT="${OS}-${OSREV}-${ARCHITECTURE}"
else
OS=SunOS
fi
fi

if test "$OS" = Linux ; then
PROCESSOR=`${UNAME} -m`
ARCHITECTURE=unknown
if test "${PROCESSOR}" = sparc -o "${PROCESSOR}" = SPARC ; then
ARCHITECTURE=sparc
BINFMT=Linux-sparc
fi
if test `echo $PROCESSOR|${CUT} -c1` = i ; then
ARCHITECTURE=i386
BINFMT=Linux-i386
fi
if test `echo $PROCESSOR` = alpha ; then
ARCHITECTURE=alpha
BINFMT=Linux-alpha
fi
if test `echo $PROCESSOR` = x86_64 ; then
ARCHITECTURE=x86_64
BINFMT=Linux-x86_64
fi
if test "${ARCHITECTURE}" = "unknown" -a "${PROJECT[]_IGNORE_ARCH}" = ""; then
AC_MSG_RESULT(OS: ${OS} / hardware: ${PROCESSOR})
AC_MSG_RESULT(Sorry - this architecture is currently not supported...)
CF_ERROR
fi
fi

if test "${OS}" = IRIX64 ; then
OS=IRIX
fi

if test "${OS}" = IRIX ; then
BINFMT=IRIX-${OSREV}
fi

if test "${OS}" = OSF1 ; then
BINFMT="OSF1-${OSREV}"
PROCESSOR=`${UNAME} -m`
ARCHITECTURE=unknown
if test `echo $PROCESSOR` = alpha ; then
ARCHITECTURE=alpha
fi
fi

if test "${OS}" = Darwin ; then
BINFMT="Darwin-${OSREV}"
PROCESSOR=`${UNAME} -p`
ARCHITECTURE=`${UNAME} -m`
fi

if test "`echo $OS | ${CUT} -d_ -f1`" = "CYGWIN" ; then
OS="CYGWIN"
PROJECT[]_NO_XDR=true
fi

if test "$OS" != Linux -a "$OS" != Solaris -a "$OS" != IRIX \
-a  "$OS" != OSF1 -a "$OS" != FreeBSD -a "$OS" != "CYGWIN" \
-a "${OS}" != Darwin -a "${PROJECT[]_IGNORE_ARCH}" = "" ; then
AC_MSG_RESULT(Sorry - your OS ($OS) is currently not supported...)
CF_ERROR
fi

dnl
dnl 	create OS defines in config.h:
dnl
if test "${OS}" = Linux ; then
AC_DEFINE(PROJECT[]_OS_LINUX,LINUX)
fi
if test "${OS}" = Solaris ; then
AC_DEFINE(PROJECT[]_OS_SOLARIS,SOLARIS)
fi
if test "${OS}" = IRIX ; then
AC_DEFINE(PROJECT[]_OS_IRIX,IRIX)
fi
if test "${OS}" = OSF1 ; then
AC_DEFINE(PROJECT[]_OS_OSF1,OSF1)
fi
if test "${OS}" = FreeBSD ; then
AC_DEFINE(PROJECT[]_OS_FREEBSD,FREEBSD)
fi
if test "${OS}" = Darwin ; then
AC_DEFINE(PROJECT[]_OS_DARWIN,DARWIN)
fi

dnl
dnl		create ARCHITECTURE defines
dnl
if test "$ARCHITECTURE" = sparc ; then
AC_DEFINE(PROJECT[]_ARCH_SPARC,SPARC)
fi
if test "$ARCHITECTURE" = i386 ; then
AC_DEFINE(PROJECT[]_ARCH_I386,I386)
fi
if test "$ARCHITECTURE" = mips ; then
AC_DEFINE(PROJECT[]_ARCH_MIPS,MIPS)
fi
if test "$ARCHITECTURE" = alpha ; then
AC_DEFINE(PROJECT[]_ARCH_ALPHA,ALPHA)
fi

AC_MSG_RESULT($OS $OSREV (BINFMT=$BINFMT))

dnl
dnl	some definitions the depend solely on the OS
dnl
SHARED_LIB_SUFFIX=so
if test "${OS}" = HP-UX ; then
SHARED_LIB_SUFFIX=sl
fi
if test "${OS}" = Darwin ; then
SHARED_LIB_SUFFIX=dylib
fi
if test "${OS}" = CYGWIN ; then
SHARED_LIB_SUFFIX=dll
fi
AC_SUBST(SHARED_LIB_SUFFIX)
])


dnl
dnl
dnl   Declare compiler search order
dnl     1) look for compiler defined in configure
dnl     2) look for vendor supplied compilers (CC)
dnl     3) check for g++, egcs, eg++, gcc
dnl   Except for Solaris, where the vendor supplied compiler
dnl   CC (at least releases 5.0 and below) is not usable.
AC_DEFUN(CF_SEARCH_CXX,[
CXX_NAME=""
case "${OS}" in
Solaris )		CXX_SEARCH_ORDER="g++ CC ";;
IRIX ) 			CXX_SEARCH_ORDER="CC g++ ";;
OSF1 )			CXX_SEARCH_ORDER="cxx CC g++ ";;
* ) 				CXX_SEARCH_ORDER="g++ CC cxx ";;
esac

dnl
dnl   Search for the C++ compiler
dnl
AC_MSG_CHECKING(searching for C++ compiler)
if test "${CXX}" != "" ; then
if test -x "${CXX}" ; then
AC_MSG_RESULT(from the command line: ${CXX})
else
AC_PATH_PROG(CXXPATH,${CXX},no)
if test "${CXXPATH}" = no ; then
AC_MSG_RESULT()
AC_MSG_RESULT(Cannot find ${CXX}. Please add it to your PATH)
AC_MSG_RESULT(or specify an absolute path in configure.)
CF_ERROR
else
CXX=${CXXPATH}
fi
fi
else
AC_MSG_RESULT()
CXXPATH=""
while test "${CXXPATH}" = "" ; do
CXX=`echo ${CXX_SEARCH_ORDER}|${CUT} -d\  -f1`
if test _`echo ${CXX} | ${TR} -d " "`_ = __ ; then
CXXPATH="END"
fi
if test "${CXXPATH}" != "END" ; then
AC_PATH_PROG(CXXPATH,${CXX},no)
if test "${CXXPATH}" = no ; then
	CXXPATH=""
	unset ac_cv_path_CXXPATH
else
	CXX=${CXXPATH}
fi
fi

CXX_SEARCH_ORDER=`echo "${CXX_SEARCH_ORDER} " |${CUT} -d\  -f2-`
done

if test "${CXXPATH}" = "end" ; then
AC_MSG_RESULT()
AC_MSG_RESULT(Could not find a C++ compiler. Please change the settings)
AC_MSG_RESULT(of your PATH environment variable (using setenv export))
AC_MSG_RESULT(or specify an absolute path in configure by setting the variable)
AC_MSG_RESULT(CXX=<pathname> or specify the compiler by passing the option)
AC_MSG_RESULT(--with-compiler=<compiler> to configure.)
CF_ERROR
fi
fi

dnl
dnl   extract the executable name of the compiler
dnl   as default compiler name (CXX_NAME is needed
dnl   to name the default directory the libraries
dnl   reside in)
dnl

if test "${CXX_PATH}" = "" ; then
if test "${CXX}" = "" ; then
CXX_NAME=unknown
else
CXX_NAME="${CXX}"
fi
else
CXX_NAME="${CXX_PATH}"
fi

while test "`echo ${CXX_NAME}|  ${GREP} /`" != "" ; do
CXX_NAME=`echo ${CXX_NAME} |  ${CUT} -d/ -f2-`
done
])



dnl
dnl		Break up the compiler version into its major and minor release numbers
dnl
AC_DEFUN(CF_DIGEST_CXX_VERSION,[
CXX_VERSION_1=`echo ${CXX_VERSION} | ${CUT} -d. -f1`
CXX_VERSION_LENGTH=`echo ${CXX_VERSION} | sed "s/[^.]//g" | wc -c`
if test "${CXX_VERSION_LENGTH}" -ge 2 ; then
CXX_VERSION_2=`echo ${CXX_VERSION} | ${CUT} -d. -f2`
fi
if test "${CXX_VERSION_LENGTH}" -ge 3 ; then
CXX_VERSION_3=`echo ${CXX_VERSION} | ${CUT} -d. -f3`
fi
if test "${CXX_VERSION_LENGTH}" -ge 4 ; then
CXX_VERSION_4=`echo ${CXX_VERSION} | ${CUT} -d. -f4`
fi
AC_DEFINE_UNQUOTED(PROJECT[]_COMPILER_VERSION_MAJOR, ${CXX_VERSION_1})
AC_DEFINE_UNQUOTED(PROJECT[]_COMPILER_VERSION_MINOR, ${CXX_VERSION_2})
AC_DEFINE_UNQUOTED(PROJECT[]_COMPILER_VERSION_MINOR_MINOR, ${CXX_VERSION_3})
])

dnl
dnl		Check whether CXX is a GNU compiler and retrieve its
dnl			version number.
dnl
AC_DEFUN(CF_IDENTIFY_GXX,[
AC_MSG_CHECKING(for GNU compiler)
cat > /tmp/$$.conftest.c << EOF
#ifdef __GNUC__
GXX:true
#else
GXX:false
#endif
EOF

IS_GXX=`${CXX} -E /tmp/$$.conftest.c 2>/dev/null | ${EGREP} GXX|${CUT} -d: -f2|${TR} -d " "`
if test "${IS_GXX}" = "true" ; then
AC_MSG_RESULT(yes)
HAS_GPLUSPLUS=true
CXX_NAME="g++"
CXX_IDENTIFIED=true

dnl 
dnl 	Define a symbol for G++.
dnl
AC_DEFINE(PROJECT[]_COMPILER_GXX, )
AC_DEFINE(PROJECT[]_COMPILER, GXX)
else
AC_MSG_RESULT(no)
HAS_GPLUSPLUS=false
fi
${RM} /tmp/$$.conftest.c

])

AC_DEFUN(CF_GXX_OPTIONS, [
AC_MSG_CHECKING(compiler version)	
VERSION_FILE=/tmp/$$.gnu_version.C
echo "__GNUC__.__GNUC_MINOR__.__GNUC_PATCHLEVEL__" > ${VERSION_FILE}
CXX_VERSION=`${CXX} -E ${VERSION_FILE} | ${GREP} -v "^#" | ${TR} -d " "`
${RM} ${VERSION_FILE}
if test `echo ${CXX_VERSION}|${CUT} -c1-4` = "egcs" ; then
IS_EGXX=true
CXX_NAME="egcs"
CXX_VERSION=`${CXX} -v 2>&1 | grep release | ${CUT} -d\( -f2 | cut -d\) -f1 | ${SED} "s/egcs-//" | ${CUT} -d" " -f1`
VERSION_OUTPUT="egcs ${CXX_VERSION}"
CXX_COMPILER_NAME="egcs"
else
IS_EGXX=false
VERSION_OUTPUT="g++ ${CXX_VERSION}"
CXX_COMPILER_NAME="g++"
fi

AC_MSG_RESULT(${VERSION_OUTPUT})
CF_DIGEST_CXX_VERSION

if test "${CXX_VERSION_1}" -lt 2 \
-o "${CXX_VERSION_1}" = 2 -a "${CXX_VERSION_2}" -lt 95 ; then
AC_MSG_RESULT()
AC_MSG_RESULT([The version of gcc you are using is not supported by PROJECT[].])
AC_MSG_RESULT([Please update to a newer version of g++ (at least 2.95.x)])
AC_MSG_RESULT([which can be obtained from])
AC_MSG_RESULT([  ftp://gcc.gnu.org/pub/gcc/releases/index.html])
AC_MSG_RESULT([or specify a different compiler using the option --with-compiler=])
CF_ERROR
fi

dnl
dnl		For g++ 3 and above we have to link against libiberty.a in order
dnl			to get the name demangling done.
dnl
if test "${CXX_VERSION_1}" -ge 3 ; then
AC_MSG_CHECKING(whether libiberty is required)
AC_MSG_RESULT(yes)
AC_MSG_CHECKING(whether libiberty is available)
SAVE_LIBS="${LIBS}"
LIBS="${LIBS} -liberty"
HAS_LIBIBERTY=false
AC_TRY_LINK([],[], HAS_LIBIBERTY=true)
if test "${HAS_LIBIBERTY}" != true ; then
LIBS="${SAVE_LIBS}"
AC_MSG_RESULT(no)
else
AC_MSG_RESULT(yes)
fi
fi


dnl
dnl		Here go the g++-specific options
dnl
CXXFLAGS="${CXXFLAGS} -pipe"
CXX_MAKEDEPEND="${CXX}"
MAKEDEP_CXX_OPTS="-M"
CXXFLAGS_D="${CXXFLAGS_D} -Wall -W -pedantic -Wno-long-long"
CXXFLAGS_DI="${CXXFLAGS_DI} -g"
dnl
dnl	Some compiler versions have problems with -O3 unter Darwin (3.3.0),
dnl so we go back to -O1.
dnl
if test "${OS}" = "Darwin" ; then
CXXFLAGS_O="${CXXFLAGS_O} -O1 -Wall -W -pedantic -Wno-long-long -Wno-long-double"		
CXXFLAGS_D="${CXXFLAGS_D} -O1"
else
CXXFLAGS_O="${CXXFLAGS_O} -O3 -Wall -W -pedantic -Wno-long-long"
fi
MAKEDEP_CXX_SUFFIX=" >.Dependencies"

dnl  We do not need the -fPIC flag for CYGWIN and Darwin,
dnl  because its code is always position independent.
dnl  A warning is emitted otherwise.
if test "${OS}" != "CYGWIN" -a "${OS}" != "Darwin" ; then
CXXFLAGS="${CXXFLAGS} -fPIC"
fi

DYNAR="${CXX}"
if test "${OS}" = "Solaris" ; then
DYNAROPTS="${DYNAROPTS} -G -fPIC -o"
else 
if test "${OS}" = "Darwin" ; then
DYNAROPTS="${DYNAROPTS} -headerpad_max_install_names -compatibility_version ${BALL_COMPATIBILITY_VERSION} -current_version ${BALL_CURRENT_VERSION} -dynamiclib -o"
ADD_DYNAROPTS_LIBBALL="-seg1addr 0xb0000000"
ADD_DYNAROPTS_LIBVIEW="-seg1addr 0x80000000"
DYNAROPTS_PYMODULE="-headerpad_max_install_names -bundle -framework Python -o"
RANLIB="ranlib -s "
else	
DYNAROPTS="${DYNAROPTS} -shared -fPIC -o"
fi
fi

if test "${IS_EGXX}" = true; then
PROJECT[]_TYPENAME=typename
else
if test "${CXX_VERSION_1}" -gt 2 -o "${CXX_VERSION_1}" -eq 2 -a "${CXX_VERSION_2}" -ge 8 ; then
PROJECT[]_TYPENAME=typename
fi
fi
])


dnl		Check for KAI C++ (KCC)
dnl   At least under linux the damned frontend won't tell
dnl   its version number, so we try to extract the word kai
dnl   from its drivers options when called in verbose mode.
dnl   Nasty - but seems to work. Anybody with a better solution
dnl   should feel free to inform me!
dnl
AC_DEFUN(CF_IDENTIFY_KAI, [
AC_MSG_CHECKING(for KAI C++ compiler)
KAI=`${CXX} -v --version 2>&1 | sed "s/.*KAI.*/__KAI__/g" |sed "s/.*kai.*/__KAI__/g" | ${EGREP} "^__KAI__$" | sed -n 1p`
if test "${KAI}" = "__KAI__" ; then
IS_KCC=true
AC_MSG_RESULT(yes)
CXX_NAME="KAI"
CXX_IDENTIFIED=true

dnl 
dnl 	Define a symbol for KAI C++.
dnl
AC_DEFINE(PROJECT[]_COMPILER_KAI, )
AC_DEFINE(PROJECT[]_COMPILER, KAI)
else
IS_KCC=false
AC_MSG_RESULT(no)
fi

])

dnl
dnl		KAI-C++-specific options
dnl
AC_DEFUN(CF_KAI_OPTIONS, [
AC_MSG_CHECKING(compiler version)
echo "int main(){}" > conftest.C
CXX_VERSION=`${CXX} -v --version conftest.C 2>&1| ${GREP} "KAI C++ " | ${CUT} -d" " -f3`
CXX_NAME="KCC"
VERSION_OUTPUT="KAI C++ ${CXX_VERSION}"
CXX_COMPILER_NAME="KCC"

AC_MSG_RESULT(${VERSION_OUTPUT})
CF_DIGEST_CXX_VERSION

dnl   KAI C++ stores a list of instantiated templates
dnl   in directories called ti_files
dnl   make clean should remove these
TEMPLATE_DIR="ti_files"
AR="${CXX}"
DYNAR="${CXX}"
AROPTS="${AROPTS} -o"
DYNAROPTS="${DYNAROPTS} -o"
CXX_MAKEDEPEND="${CXX}"
MAKEDEP_CXX_OPTS="-M"
MAKEDEP_CXX_SUFFIX=" >.Dependencies"

dnl
dnl   Someone at KAI seems to have the need
dnl   to torture developers by introducing
dnl   a new flag for position independent code
dnl   on EVERY platform...
dnl
CXXFLAGS="${CXXFLAGS} --one_per"
if test "${OS}" = Linux ; then
CXXFLAGS="${CXXFLAGS} -fPIC"
fi
if test "${OS}" = Solaris ; then
CXXFLAGS="${CXXFLAGS} -KPIC"
fi
if test "${OS}" = IRIX ; then
CXXFLAGS="${CXXFLAGS} -KPIC"
fi

dnl   optimze as on highest level: this compiler
dnl   does a good job optimizing!
CXXFLAGS_O="${CXXFLAGS_O} +K3"

dnl   avoid high level optimization to
dnl   get debuggable code...
CXXFLAGS_D="${CXXFLAGS_D} +K0"
CXXFLAGS_DI="${CXXFLAGS_DI}"

dnl
dnl if we are running under Solaris/SPARC,
dnl KAI can produce 32 or 64 bit code
dnl
if test "${OS}" = "Solaris" -a "${ARCHITECTURE}" = sparc ; then
if test "${BINFMT_64_BIT}" = true ; then
CXX_NAME="${CXX_NAME}_V9"
LDFLAGS="${LDFLAGS} -xarch=v9"
CXXFLAGS="${CXXFLAGS} -xarch=v9"
AROPTS="${AROPTS} -xarch=v9"
DYNAROPTS="-xarch=v9 ${DYNAROPTS}"
else
CXX_NAME="${CXX_NAME}_V8"
fi
fi
])

dnl
dnl   Check for Intel C++ (icc)
dnl
AC_DEFUN(CF_IDENTIFY_INTEL, [
AC_MSG_CHECKING(for Intel C++ compiler)
ICC=`${CXX} -V 2>&1 | ${SED} -n 1p |${SED} "s/Intel(R) C Compiler.*/__INTELCC__/g"| ${SED} "s/Intel(R) C++ Compiler.*/__INTELCC__/g" | ${EGREP} "^__INTELCC__$" | sed -n 1p`
if test "${ICC}" = "__INTELCC__" ; then
IS_INTELCC=true
AC_MSG_RESULT(yes)
CXX_NAME="Intel"
CXX_IDENTIFIED=true


dnl 
dnl 	Define a symbol for Intel C++.
dnl
AC_DEFINE(PROJECT[]_COMPILER_INTEL, )
AC_DEFINE(PROJECT[]_COMPILER, INTEL)
else
IS_INTELCC=false
AC_MSG_RESULT(no)
fi
])

dnl
dnl		Set the Intel icc-specific options.
dnl
AC_DEFUN(CF_INTEL_OPTIONS,[
AC_MSG_CHECKING(compiler version)
echo "int main(){}" > conftest.C
CXX_VERSION=`${CXX} -V conftest.C 2>&1| ${GREP} "Intel(R) C++" | ${SED} -n 1p | ${SED} "s/.*Version //" | ${CUT} -d" " -f1`
if test "${CXX_VERSION}" = "" ; then
  CXX_VERSION=`${CXX} -v 2>&1| ${SED} "s/.*ersion //" | ${CUT} -d " " -f1`
fi
CXX_NAME="icc"
VERSION_OUTPUT="Intel C++ Compiler ${CXX_VERSION}"
CXX_COMPILER_NAME="icc"

AC_MSG_RESULT(${VERSION_OUTPUT})
CF_DIGEST_CXX_VERSION


DYNAROPTS="${DYNAROPTS} -shared -o"
TEMPLATE_DIR=""
AR="${CXX}"
DYNAR="${CXX}"
AROPTS="${AROPTS} -o"
CXX_MAKEDEPEND="${CXX}"
MAKEDEP_CXX_OPTS="-M"
MAKEDEP_CXX_SUFFIX=" >.Dependencies"
CXXFLAGS="${CXXFLAGS} -KPIC"

dnl
if test "${CXX_VERSION_1}" != "" -a "${CXX_VERSION_2}" != "" ; then
  if test ${CXX_VERSION_1} = 8 -a ${CXX_VERSION_2} -ge 1 -o ${CXX_VERSION_1} -ge 9 ; then
    if test `basename ${CXX}` = "icc" ; then
      AC_MSG_RESULT([WARNING: Starting with version 8.1, icpc should be used instead of icc.])
      AC_MSG_RESULT([Otherwise, linking errors will occur! Please call configure again])
      AC_MSG_RESULT([with icpc as the compiler.])
      AC_MSG_RESULT()
      CF_ERROR(Aborted.)
    fi
  fi 
  dnl
  dnl  icpc 9.0 has a problem: linking without -O0 invokes IPO, which
  dnl  doesn't seem to work ("Cannot find "("..." as error message).
  dnl  So we just turn it off.
  dnl
  if test ${CXX_VERSION_1} -ge 9 ; then
    LDFLAGS="$LDFLAGS -O0"
    DYNAROPTS="-O0 $DYNAROPTS"
  fi
fi



dnl   optimze as on highest level: this compiler
CXXFLAGS_O="${CXXFLAGS_O} -O1"

dnl   avoid high level optimization to
dnl   get debuggable code...
CXXFLAGS_D="${CXXFLAGS_D} -O0 -g -w1"
CXXFLAGS_DI="${CXXFLAGS_DI}"
])

dnl
dnl   check for the Digital C++ compiler
dnl
dnl
AC_DEFUN(CF_IDENTIFY_COMPAQ,[
AC_MSG_CHECKING(for Digital/Compaq C++ compiler)
DIGITAL_CXX=`${CXX} -V 2>/dev/null | ${GREP} "C++" | ${CUT} -d" " -f-2`
if test "${DIGITAL_CXX}" = "DIGITAL C++" -o "${DIGITAL_CXX}" = "Compaq C++"; then
IS_DIGITALCXX=true
AC_MSG_RESULT(yes)
CXX_NAME="Compaq"
CXX_IDENTIFIED=true

dnl 
dnl 	Define a symbol for Compaq C++.
dnl
AC_DEFINE(PROJECT[]_COMPILER_COMPAQ, )
AC_DEFINE(PROJECT[]_COMPILER, COMPAQ)
else
IS_DIGITALCXX=false
AC_MSG_RESULT(no)
fi

])

AC_DEFUN(CF_COMPAQ_OPTIONS, [
AC_MSG_CHECKING(compiler version)
echo "int main(){}" > conftest.C
CXX_VERSION=`${CXX} -V  2>/dev/null| ${GREP} "C++" | ${SED} "s/^.*C++ //" | ${CUT} -d" " -f1 | ${TR} -d V | ${TR} "-" "."`
CXX_NAME=`${CXX} -V | ${GREP} "C++" | ${CUT} -d" " -f1`
VERSION_OUTPUT="${CXX_NAME} C++ ${CXX_VERSION}"
CXX_COMPILER_NAME="Digital"

AC_MSG_RESULT(${VERSION_OUTPUT})
CF_DIGEST_CXX_VERSION

if test "${CXX_VERSION_1}" -lt 6 -o "${CXX_VERSION_1}" -eq 6 -a "${CXX_VERSION_2}" -lt 2 ; then
AC_MSG_RESULT()
AC_MSG_RESULT(Your version of Digital/Compaq C++ does not provide all)
AC_MSG_RESULT(ANSI C++ features required by PROJECT[].)
AC_MSG_RESULT(Please upgrade to release 6.2 or above.)
CF_ERROR
fi

TEMPLATE_DIR="cxx_rep"
AR="ar"
DYNAR="${CXX}"
AROPTS="${AROPTS} -o"
DYNAROPTS="${DYNAROPTS} -shared -nocxxstd -ptr \$(PROJECT[]_PATH)/source/cxx_rep -o"
CXX_MAKEDEPEND="${CXX}"
MAKEDEP_CXX_OPTS="-M -noimplicit_include"
MAKEDEP_CXX_SUFFIX=" >.Dependencies"

dnl 
dnl  CXX 6.2 does not provide the -nopure_cname flag
if test "${CXX_VERSION_2}" -lt 3 ; then
CXXFLAGS="${CXXFLAGS} -ieee"
else
if test "${CXX_VERSION_2}" -ge 5 ; then
dnl
dnl  Starting with cxx 6.5, we had some trouble with
dnl  floating point accuracy -- that should take care of it.
CXXFLAGS="${CXXFLAGS} -ieee -nopure_cname"
else
CXXFLAGS="${CXXFLAGS} -ieee -nopure_cname"
fi
fi

LIB_CXXFLAGS="${LIB_CXXFLAGS} -ptr \$(PROJECT[]_PATH)/source/cxx_rep"
CXXFLAGS_O="${CXXFLAGS_O} -O3"

CXXFLAGS_D="${CXXFLAGS_D}"
CXXFLAGS_DI="${CXXFLAGS_DI} -g"

dnl   Problem with linux headers:
dnl   cannot use -std strict_ansi since the socket headers
dnl   cause an error #280
if test "${OS}" != "Linux" ; then
CXXFLAGS="${CXXFLAGS} -std strict_ansi"
MAKEDEP_CXX_OPTS="${MAKEDEP_CXX_OPTS} -std strict_ansi"
fi
])

AC_DEFUN(CF_IDENTIFY_SGI, [
AC_MSG_CHECKING(for SGI/MipsPro C++ compiler)
SGI_CXX=`${CXX} -version -n32 2>&1 | ${GREP} "MIPSpro" | ${CUT} -d" " -f1`
if test "${SGI_CXX}" = "MIPSpro"; then
IS_MIPSPRO=true
AC_MSG_RESULT(yes)
CXX_NAME="MIPSpro"
CXX_IDENTIFIED=true

dnl 
dnl 	Define a symbol for SGI C++.
dnl
AC_DEFINE(PROJECT[]_COMPILER_MIPSPRO, )
AC_DEFINE(PROJECT[]_COMPILER, MIPSPRO)
else
IS_MIPSPRO=false
AC_MSG_RESULT(no)
fi
])

AC_DEFUN(CF_MIPSPRO_OPTIONS, [
AC_MSG_CHECKING(compiler version)
CXX_VERSION_STRING=`${CXX} -n32 -version 2>&1 | ${EGREP} ersion`
if test "${CXX_VERSION_STRING}" = "" ; then
CF_BASENAME(${CXX})
CXX_VERSION="${TMP__NAME}"
CXX_VERSION_OUTPUT="${CXX_VERSION} (unknown version)"
else
CXX_VERSION=`echo ${CXX_VERSION_STRING} | ${SED} "s/^.*ersion //g"`
CXX_COMPILER_NAME=`echo ${CXX_VERSION_STRING} | ${CUT} -d\  -f1`
CXX_VERSION_OUTPUT="${CXX_VERSION} (${CXX_COMPILER_NAME})"
fi
CF_DIGEST_CXX_VERSION	
AC_MSG_RESULT(${CXX_VERSION_OUTPUT})

dnl  set the name for the template repository
dnl
TEMPLATE_DIR="ii_files"

dnl  set the default binary format (if none selected)
dnl
if test "${BINFMT_64_BIT}" = true ; then
IRIX_BINFMT=64
CXX_NAME="${CXX_NAME}_64"
else
IRIX_BINFMT=N32
CXX_NAME="${CXX_NAME}_N32"
fi

PROJECT[]_TYPENAME=typename

dnl
dnl     a version above 7.2 is required
dnl
if test "${CXX_VERSION_1}" -lt 7\
		-o "${CXX_VERSION_1}" -eq 7 -a "${CXX_VERSION_2}" = 10\
		-o "${CXX_VERSION_1}" -eq 7 -a "${CXX_VERSION_2}" = 20\
		-o "${CXX_VERSION_1}" -eq 7 -a "${CXX_VERSION_2}" -lt 2; then
AC_MSG_RESULT()
AC_MSG_RESULT(MipsPro CC version 7.30 or above is required. Please update your compiler.)
CF_ERROR
fi

AR=${CXX}
AROPTS="${AROPTS} -ar -o"
DYNAR=${CXX}
DYNAROPTS="${DYNAROPTS} -LANG:std -shared -quickstart_info -no_unresolved -o"

dnl  issue a warning about an old compiler with a broken ostream implementation: reopening a fstream
dnl  and writing to it will omit the first 16k written to the stream. Nasty, but confirmed with SGI
dnl  and fixed in 7.3.1.1m
if test "${CXX_VERSION_1}" = 7 -a "${CXX_VERSION_2}" = 3 \
		-a "${CXX_VERSION_3}" = 1m -a "${CXX_VERSION_4}" = ""\
		-o "${CXX_VERSION_1}" = 7 -a "{CXX_VERSION_2}" = 30 ; then
COMMENTS="${COMMENTS}\nPlease take care - this version of SGI CC (7.3.1m) contains serious bugs\n"
COMMENTS="${COMMENTS}in its implementaion of fstream/iostream. This may lead to strange behaviour\n"
COMMENTS="${COMMENTS}and causes PDBFile_test to fail. Please update your compiler.\n\n"
AC_MSG_RESULT(${COMMENTS})
ADDITIONAL_COMMENTS="${ADDITIONAL_COMMENTS}${COMMENTS}"
COMMENTS=""
fi

CXX_MAKEDEPEND="${CXX}"
MAKEDEP_CXX_OPTS="-M 2>/dev/null"
MAKEDEP_CXX_SUFFIX=" >.Dependencies"

if test "${IRIX_BINFMT}" = 64 ; then
DEF_BOOL=false
CXXFLAGS="$CXXFLAGS -64 -LANG:std"
DYNAROPTS="-64 ${DYNAROPTS}"
CXXFLAGS_O="${CXXFLAGS_O} -O3 -OPT:Olimit=60000 -multigot -G 5 -DEBUG:woff=3333"
CXXFLAGS_D="${CXXFLAGS_D} -fullwarn -multigot -G 5 -DEBUG:woff=1375,3201,1424,3333,1110,1209"
CXXFLAGS_DI="${CXXFLAGS_DI} -g"
LDFLAGS="$LDFLAGS -64 -LANG:std"
AC_DEFINE(IRIX64,)
fi
if test "${IRIX_BINFMT}" = N32 ; then
DEF_BOOL=false
CXXFLAGS="$CXXFLAGS -n32 -LANG:std"
DYNAROPTS="-n32 ${DYNAROPTS}"
CXXFLAGS_O="${CXXFLAGS_O} -O3 -OPT:Olimit=60000 -multigot -G 5 -DEBUG:woff=3333"
CXXFLAGS_D="${CXXFLAGS_D} -fullwarn -multigot -G 5 -DEBUG:woff=1375,3201,1424,3333,1110,1209"
CXXFLAGS_DI="${CXXFLAGS_DI} -g"
LDFLAGS="$LDFLAGS -n32 -LANG:std"
AC_DEFINE(IRIX32,)
fi

dnl
dnl  -O3 requires -IPA for linking
dnl
if test "${DEBUG}" = yes ; then
DYNAROPTS="-IPA ${DYNAROPTS}"
fi

AC_DEFINE(MIPS,)
AC_DEFINE(IRIX,)
])

AC_DEFUN(CF_IDENTIFY_SUN, [
AC_MSG_CHECKING(for SUN WorkShop/Forte C++ compiler)
SUN_CXX=`${CXX} -V 2>&1 | ${EGREP} -e "(Sun)|(WorkShop)|(Forte)"`
if test "${SUN_CXX}" != ""; then
IS_SUNCC=true
AC_MSG_RESULT(yes)
CXX_NAME="SunCC"
CXX_IDENTIFIED=true

dnl 
dnl 	Define a symbol for SUNPro C++.
dnl
AC_DEFINE(PROJECT[]_COMPILER_SUNPRO)
AC_DEFINE(PROJECT[]_COMPILER, SUNPRO)
else
IS_SUNCC=false
AC_MSG_RESULT(no)
fi
])

AC_DEFUN(CF_SUNCC_OPTIONS, [

AC_MSG_CHECKING(compiler version for Sun C++)
CXX_VERSION=`echo ${SUN_CXX} | ${SED} "s/^.*C++ //" | ${CUT} -d" " -f1`
changequote(<<,>>)
CXX_VERSION_TEST=`echo ${CXX_VERSION} | ${SED} s/\^\\[0-9\\.\\]*[a-zA-Z\\.]*//g`
changequote([,])
if test "${CXX_VERSION_TEST}" != "" ; then
CF_BASENAME(${CXX})
CXX_VERSION="${TMP__NAME}"
CXX_VERSION_OUTPUT="${CXX_VERSION} (unknown version)"
else
CXX_VERSION_OUTPUT="${CXX_VERSION}"
CF_DIGEST_CXX_VERSION	
fi
AC_MSG_RESULT(${CXX_VERSION_OUTPUT})

dnl
dnl   Make sure we use at least Workshop 6U2 (C++ 5.3)
dnl   (SUNPro < 6 is a mess - hasn't even heard of ANSI C++!)
dnl
if test "${CXX_VERSION_1}" -lt 5 ; then
AC_MSG_RESULT()
AC_MSG_RESULT(PROJECT[] requires an ANSI C++ compliant compiler)
AC_MSG_RESULT(SUNPro compilers are (mostly) ANSI compliant for version 5.3 and above)
AC_MSG_RESULT(Please upgrade your compiler!)
CF_ERROR
fi

AC_DEFINE(SOLARIS,)

TEMPLATE_DIR="SunWS_cache"

dnl  a nasty bug in SUNPro CC 5.3 causes trouble
dnl  with the function templates in amberNonBonded.C
AC_DEFINE(PROJECT[]_MUST_CAST_TEMPLATE_FUNCTION_ARGS,)

dnl  set the default binary format (if none selected)
dnl
if test "${BINFMT_64_BIT}" = true ; then
SUN_BINFMT=V9
CXX_NAME="${CXX_NAME}_V9"
LDFLAGS="${LDFLAGS} -xarch=v9"
CXXFLAGS="${CXXFLAGS} -xarch=v9"
AROPTS="${AROPTS} -xarch=v9"
DYNAROPTS="-xarch=v9 ${DYNAROPTS}"
else
SUN_BINFMT=V8
CXX_NAME="${CXX_NAME}_V8"
fi

DEF_BOOL=true
AR="${CXX}"
DYNAR="${CXX}"
AROPTS="${AROPTS} -xar -KPIC -o"
DYNAROPTS="${DYNAROPTS} -pto -G -KPIC -o"
NONLIB_CXXFLAGS="-pto"
LIB_CXXFLAGS=""
CXX_MAKEDEPEND="${CXX}"
MAKEDEP_CXX_OPTS="-xM1"
MAKEDEP_CXX_SUFFIX=" >.Dependencies"

AC_DEFINE(PROJECT[]_NO_INLINE_FUNCTIONS,)

CXXFLAGS="${CXXFLAGS} -KPIC"
CXXFLAGS_O="${CXXFLAGS_O} -xO5"
CXXFLAGS_D="${CXXFLAGS_D}"
CXXFLAGS_DI="${CXXFLAGS_DI} -g"
])

dnl
dnl   Assemble the complete compiler name by adding
dnl   the release numbers (if known) of the compiler
dnl
AC_DEFUN(CF_BUILD_FULL_CXX_NAME, [
AC_MSG_CHECKING(standardized compiler name)
if test "${CXX_VERSION_1}" != "" ; then
CXX_NAME="${CXX_NAME}_${CXX_VERSION_1}"
if test "${CXX_VERSION_2}" != "" ; then
CXX_NAME="${CXX_NAME}.${CXX_VERSION_2}"
if test "${CXX_VERSION_3}" != "" ; then
CXX_NAME="${CXX_NAME}.${CXX_VERSION_3}"
if test "${CXX_VERSION_4}" != "" ; then
	CXX_NAME="${CXX_NAME}.${CXX_VERSION_4}"
fi
fi
fi
fi

AC_MSG_RESULT(${CXX_NAME})
])

dnl
dnl
dnl   checking for DEBUG-Flag
dnl
AC_DEFUN(CF_CHECK_DEBUG_FLAG, [
AC_MSG_CHECKING(for DEBUG flag)
if test "$DEBUG" != "" ; then
dnl   define a debug flag and prevent the compilation of
dnl   inline functions by defining PROJECT[]_NO_INLINE_FUNCTIONS
dnl   (see COMMON/debug.h)
if test "$DEBUG" = true ; then
dnl  if debug information is also required, add the corresponding flag
dnl
if test "${DEBUG_INFO}" = true -a "$CXXFLAGS_DI" != "" ; then
CXXFLAGS_D="${CXXFLAGS_D} ${CXXFLAGS_DI}"
fi
AC_DEFINE(PROJECT[]_DEBUG,)
AC_DEFINE(PROJECT[]_NO_INLINE_FUNCTIONS,)
AC_MSG_RESULT(enabled)
CPP_MODE_FLAGS="${CXXFLAGS_D}"
CPP_MODE_FLAGS_NO_OPTIMIZATION="${CXXFLAGS_D}"
else
AC_MSG_RESULT(disabled)
CPP_MODE_FLAGS="${CXXFLAGS_O}"
CPP_MODE_FLAGS_NO_OPTIMIZATION=""
fi
else
AC_MSG_RESULT(disabled)
CPP_MODE_FLAGS="${CXXFLAGS_O}"
CPP_MODE_FLAGS_NO_OPTIMIZATION=""
fi
])

dnl
dnl   check for endianness of the architecture
dnl
dnl
AC_DEFUN(CF_C_BIGENDIAN, [
AC_MSG_CHECKING(for byte order)
AC_TRY_RUN(
[
#include <iostream>
#include <fstream>
int main(int, char**)
{
] ${PROJECT[]_SIZE_TYPE} endian_one = 1; [
std::ofstream os("config.endian.log", std::ios::out);


if (*(char*)&endian_one == '\001')
{
// big endian
os << "LITTLE" << std::endl;
}
else
{
// little endian
os << "BIG" << std::endl;
}
os.close();

return 0;
}
],
PROJECT[]_ENDIAN_TEST=true,
DUMMY=0,
DUMMY=0
)
if test "${PROJECT[]_ENDIAN_TEST+set}" != set ; then
AC_MSG_RESULT(<cannot determine>)
CF_ERROR
else
dnl
dnl read the result of the endian test from the file
dnl and delete the file
dnl
ENDIAN_TYPE=`${CAT} config.endian.log`
${RM} config.endian.log 2>/dev/null
if test "${ENDIAN_TYPE}" = "LITTLE" ; then
PROJECT[]_LITTLE_ENDIAN=true
AC_DEFINE(PROJECT[]_LITTLE_ENDIAN, true)
AC_MSG_RESULT(little endian)
else
if test "${ENDIAN_TYPE}" = "BIG" ; then
PROJECT[]_BIG_ENDIAN=true
AC_DEFINE(PROJECT[]_BIG_ENDIAN, true)
AC_MSG_RESULT(big endian)
else
AC_MSG_RESULT(<cannot determine>)
CF_ERROR
fi
fi
fi
])


dnl
dnl   check for limits header (class numeric limits, to be precise)
dnl
dnl
AC_DEFUN(CF_CHECK_NUM_LIMITS, [
AC_MSG_CHECKING(for numeric_limits class)
AC_TRY_COMPILE(
[
#include <limits>
],
[
float f = std::numeric_limits<float>::min();
],
PROJECT[]_HAS_NUMERIC_LIMITS=true
)
if test "${HAS_NUMERIC_LIMITS}" = true ; then
AC_MSG_RESULT(available)
AC_DEFINE(PROJECT[]_HAS_NUMERIC_LIMITS)
else
AC_MSG_RESULT(not available)
PROJECT[]_HAS_NUMERIC_LIMITS=false

dnl
dnl  we didn't find a numeric limits class, so we implement
dnl  it by ourselves. Try to figure out whether we only need
dnl  limts.h or also float.h (sometimes FLT_MIN seems to be missing)
dnl
AC_MSG_CHECKING(whether float.h has to be included)
AC_TRY_COMPILE(
[
#include <limits.h>
],
[
float a = FLT_MAX;
float b = FLT_MIN;
float c = DBL_MAX;
float d = DBL_MIN;
],
PROJECT[]_HAS_FLOAT_H=false
)
if test "${PROJECT[]_HAS_FLOAT_H}" != false ; then
AC_TRY_COMPILE(
[
	#include <float.h>
],
[
	float a = FLT_MAX;
	float b = FLT_MIN;
	float c = DBL_MAX;
	float d = DBL_MIN;
],
PROJECT[]_HAS_FLOAT_H=true
)
fi
if test "${PROJECT[]_HAS_FLOAT_H+set}" != set ; then
AC_MSG_RESULT()
AC_MSG_RESULT(limits.h seems to be corrupt or float.h is missing!)
AC_MSG_RESULT()
else
if test "${PROJECT[]_HAS_FLOAT_H}" = true ; then
AC_MSG_RESULT(yes)
AC_DEFINE(PROJECT[]_HAS_FLOAT_H)
else
AC_MSG_RESULT(no)
fi
fi
fi
${RM} /tmp/$$.conftest.c a.out 2>/dev/null
])


dnl
dnl   check for template arguments needed for friends of template
dnl   classses. Some compilers require "<>" after the method name,
dnl   others don't - so let's find it out!
dnl
AC_DEFUN(CF_CHECK_TPL_NULL_ARGS, [
AC_MSG_CHECKING(for null template arguments)
PROJECT[]_NULL_TEMPLATE_ARGS="NULL"
AC_TRY_COMPILE(
[
template <typename T>
class A
{
public:
friend bool operator == <> (const A&, const A&);
};
],
[
],
PROJECT[]_NULL_TEMPLATE_ARGS="<>")
if test "${PROJECT[]_NULL_TEMPLATE_ARGS}" = "NULL" ; then
AC_TRY_COMPILE(
[
template <typename T>
class A
{
	public:
	friend bool operator == (const A&, const A&);
};
],
[
],
PROJECT[]_NULL_TEMPLATE_ARGS="")
fi
AC_MSG_RESULT(\"$PROJECT[]_NULL_TEMPLATE_ARGS\")
if test "${PROJECT[]_NULL_TEMPLATE_ARGS}" = "NULL" ; then
AC_MSG_RESULT(could not find a suitable argument for null templates)
CF_ERROR
fi
])

dnl
dnl		Check whether the compiler allows parameterization oftemplate functions
dnl		with inline functions (SGI CC has a problem with that)
dnl
AC_DEFUN(CF_CHECK_INLINE_TPL_ARGS, [
AC_MSG_CHECKING(for inline template function arguments)
PROJECT[]_HAS_INLINE_TPL_ARGS=no
AC_TRY_COMPILE(
[
template <int i>
inline double foo(double x){ return i * x; }

typedef double (*Function)(double);

template <Function F>
inline double bar(double x) { return F(x); }
],
[
double d = bar< foo<3> >(2.0);
],
PROJECT[]_HAS_INLINE_TPL_ARGS=true
)
if test "${PROJECT[]_HAS_INLINE_TPL_ARGS}" = true ; then
AC_MSG_RESULT(yes)
AC_DEFINE(PROJECT[]_HAS_INLINE_TPL_ARGS)
else
AC_MSG_RESULT(no)
fi
])

dnl
dnl   check for ANSI compliant <iostream>
dnl   We need this for the base classes (ios vs. basic_ios<char>) in socket.h/C
dnl
AC_DEFUN(CF_CHECK_ANSI_IOSTREAM, [
AC_MSG_CHECKING(for ANSI compliant iostream)
PROJECT[]_HAS_ANSI_IOSTREAM=no
AC_TRY_COMPILE(
[
#include <iostream>
class A:public std::iostream
{
A():std::basic_ios<char>(0),std::iostream(0)
{}
};
],
[
],
PROJECT[]_HAS_ANSI_IOSTREAM=yes
)
AC_MSG_RESULT($PROJECT[]_HAS_ANSI_IOSTREAM)
])

dnl
dnl		check if we can use <sstream> or if we need to support old
dnl		style strstream
dnl
AC_DEFUN(CF_CHECK_HAS_SSTREAM, [
AC_MSG_CHECKING(for sstream headers)
PROJECT[]_HAS_SSTREAM=no
AC_TRY_COMPILE(
[
#include <sstream>
class A : public std::stringbuf
{
A() : std::stringbuf()
{}
};
],
[
],
PROJECT[]_HAS_SSTREAM=yes
)
AC_MSG_RESULT($PROJECT[]_HAS_SSTREAM)
])

dnl
dnl   check for ANSI or ARM style access modification
dnl   either (ARM style) Base::foo or (ANSI style) using Base::foo
dnl
AC_DEFUN(CF_CHECK_ARM_ACCESS_MODIFICATION, [
AC_MSG_CHECKING(for ANSI or ARM style access modification)
PROJECT[]_CFG_USING_METHOD_DIRECTIVE=none
AC_TRY_COMPILE(
[
class A
{
protected: void foo(){};
};

class B : public A
{
public: using A::foo;
};
],
[
B b;
b.foo();
],
PROJECT[]_CFG_USING_METHOD_DIRECTIVE=ANSI
)
if test ${PROJECT[]_CFG_USING_METHOD_DIRECTIVE} = none ; then
AC_TRY_COMPILE(
[
class A
{
	protected: void foo(){};
};

class B : public A
{
	public: A::foo;
};
],
[
B b;
b.foo();
],
PROJECT[]_CFG_USING_METHOD_DIRECTIVE=ARM
)
fi
AC_MSG_RESULT(${PROJECT[]_CFG_USING_METHOD_DIRECTIVE})
if test ${PROJECT[]_CFG_USING_METHOD_DIRECTIVE} = ANSI ; then
AC_DEFINE(PROJECT[]_CFG_USING_METHOD_DIRECTIVE)
fi
if test ${PROJECT[]_CFG_USING_METHOD_DIRECTIVE} = none ; then
AC_MSG_RESULT()
AC_MSG_RESULT([Compiler does not understand ARM or ANSI style method access modification.])
AC_MSG_RESULT([Please specify a different compiler (e.g. g++ 2.95.2) using the option])
AC_MSG_RESULT([--with-compiler=<compiler>.])
CF_ERROR
fi
])

dnl
dnl	Check for the sizes of the tzpes and define a set of portable
dnl		types of different word lengths.
dnl	This might profit very nicely from stdint.h at some point.
dnl
AC_DEFUN(CF_GET_TYPE_SIZES, [
dnl
dnl   check for the size of int and pointers (may cause trouble on 64 bit architectures)
dnl   we define the type PointerInt (in COMMON/global.h) according to the macro
dnl   PROJECT[]_POINTERSIZE_UINT (which is set here)
dnl   We also define a 64 bit unsigned numeric type. All pointers that are read or written
dnl   in persistence-related methods use this type to ensure compatibility between 32 and
dnl   64bit PROJECT[] versions.
dnl   missing: usage of the result of AC_TYPE_SIZE_T
dnl
AC_CHECK_SIZEOF(char, 4)
AC_CHECK_SIZEOF(int, 4)
AC_CHECK_SIZEOF(long, 4)
AC_CHECK_SIZEOF(long long, 8)
AC_CHECK_SIZEOF(size_t, 4)
AC_CHECK_SIZEOF(void*, 4)
SIZEOF_INT=$ac_cv_sizeof_int
SIZEOF_CHAR=$ac_cv_sizeof_char
SIZEOF_LONG=$ac_cv_sizeof_long
SIZEOF_SIZE_T=$ac_cv_sizeof_size_t
SIZEOF_VOID_P=$ac_cv_sizeof_voidp
SIZEOF_UINT=$ac_cv_sizeof_int
SIZEOF_ULONG=$ac_cv_sizeof_long
SIZEOF_ULONGLONG=$ac_cv_sizeof_long_long

AC_DEFINE_UNQUOTED(PROJECT[]_CHAR_SIZE, ${SIZEOF_CHAR})
AC_DEFINE_UNQUOTED(PROJECT[]_INT_SIZE, ${SIZEOF_INT})
AC_DEFINE_UNQUOTED(PROJECT[]_LONG_SIZE, ${SIZEOF_LONG})
AC_DEFINE_UNQUOTED(PROJECT[]_SIZE_T_SIZE, ${SIZEOF_SIZE_T})
AC_DEFINE_UNQUOTED(PROJECT[]_POINTER_SIZE, ${SIZEOF_VOID_P})
AC_DEFINE_UNQUOTED(PROJECT[]_UINT_SIZE, ${SIZEOF_UINT})
AC_DEFINE_UNQUOTED(PROJECT[]_ULONG_SIZE, ${SIZEOF_ULONG})
AC_DEFINE_UNQUOTED(PROJECT[]_ULONGLONG_SIZE, ${SIZEOF_ULONGLONG})
dnl
dnl  define an unsigned type that can hold 64 bit pointers
dnl
if test "${SIZEOF_UINT}" = 8; then
PROJECT[]_64BIT_UINT="unsigned int"
else
if test "${SIZEOF_ULONG}" = 8; then
PROJECT[]_64BIT_UINT="unsigned long"
else
if test "${SIZEOF_ULONGLONG}" = 8 ; then
PROJECT[]_64BIT_UINT="unsigned long long"
else
AC_MSG_RESULT()
AC_MSG_RESULT(cannot find appropriate numeric type for 64bit unsigned int)
CF_ERROR
fi
fi
fi
AC_DEFINE_UNQUOTED(PROJECT[]_64BIT_UINT_TYPE, ${PROJECT[]_64BIT_UINT})

dnl
dnl define a 32 bit type for Size and Index
dnl
if test "${SIZEOF_UINT}" = "${SIZEOF_VOID_P}" ; then
PROJECT[]_POINTER_TYPE="unsigned int"
else
if test "${SIZEOF_ULONG}" = "${SIZEOF_VOID_P}" ; then
PROJECT[]_POINTER_TYPE="unsigned long"
else
AC_MSG_RESULT()
AC_MSG_RESULT(cannot find appropriate integer type of same size as void*)
CF_ERROR
fi
fi
AC_DEFINE_UNQUOTED(PROJECT[]_POINTERSIZEUINT_TYPE, ${PROJECT[]_POINTER_TYPE})

dnl
dnl define a (true) pointer size int for several conversion issues
dnl (e.g. PersistenceManager) - just a define, no real typedef
dnl since for internal use only!
dnl
if test "${SIZEOF_INT}" = 4 ; then
PROJECT[]_INDEX_TYPE="int"
PROJECT[]_SIZE_TYPE="unsigned int"
else
if test "${SIZEOF_LONG}" = 4 ; then
PROJECT[]_INDEX_TYPE="long"
PROJECT[]_SIZE_TYPE="unsigned long"
else
AC_MSG_RESULT()
AC_MSG_RESULT(cannot find appropriate numeric type for 32bit int)
CF_ERROR
fi
fi
AC_DEFINE_UNQUOTED(PROJECT[]_SIZE_TYPE, ${PROJECT[]_SIZE_TYPE})
AC_DEFINE_UNQUOTED(PROJECT[]_INDEX_TYPE, ${PROJECT[]_INDEX_TYPE})

dnl  define 64 bit signed/unsigned type
if test "${SIZEOF_ULONG}" = "8" ; then
PROJECT[]_ULONG64_TYPE="unsigned long"
PROJECT[]_LONG64_TYPE="long"
AC_DEFINE(PROJECT[]_64BIT_ARCHITECTURE)
else
if test "${SIZEOF_ULONGLONG}" = "8" ; then
PROJECT[]_ULONG64_TYPE="unsigned long long"
PROJECT[]_LONG64_TYPE="long long"			
else
AC_MSG_RESULT()
AC_MSG_RESULT(cannot find unsigned 64bit type.)
CF_ERROR
fi
fi
AC_DEFINE_UNQUOTED(PROJECT[]_ULONG64_TYPE, ${PROJECT[]_ULONG64_TYPE})
AC_DEFINE_UNQUOTED(PROJECT[]_LONG64_TYPE, ${PROJECT[]_LONG64_TYPE})

dnl
dnl Check for size of Complex type
dnl
PROJECT[]_COMPLEX_PRECISION=float
AC_MSG_CHECKING(for Complex type precision)
if test "${enable_double_cplx}" = yes ; then
PROJECT[]_COMPLEX_PRECISION=double
fi
AC_MSG_RESULT(${PROJECT[]_COMPLEX_PRECISION})
])


dnl
dnl   check whether the <regex.h> header exists
dnl
AC_DEFUN(CF_CHECK_REGEX_H, [
AC_CHECK_HEADER(regex.h, HAS_REGEX_H=true, HAS_REGEX_H=false)
if test "${HAS_REGEX_H}" = "false" ; then
AC_CHECK_HEADER(regexp.h, HAS_REGEX_H=true, HAS_REGEX_H=false)
AC_DEFINE(PROJECT[]_HAS_REGEXP_H,)
else
AC_DEFINE(PROJECT[]_HAS_REGEX_H,)
fi
if test "${HAS_REGEX_H}" = "false" ; then
AC_MSG_RESULT()
AC_MSG_RESULT([Regular expression headers regex.h not found!])
AC_MSG_RESULT([If you do not have this header on your system,])
AC_MSG_RESULT([please install the GNU regexp package from])
AC_MSG_RESULT()
AC_MSG_RESULT([  ftp://ftp.gnu.org/gnu/regex/regex-0.12.tar.gz])
CF_ERROR
fi
])

dnl
dnl   Check whether ieeefp.h does really exist.
dnl
AC_DEFUN(CF_CHECK_IEEEFP_H, [
AC_CHECK_HEADERS(ieeefp.h,
[PROJECT[]_HAS_IEEEFP_H=true],
[PROJECT[]_HAS_IEEEFP_H=false])
if test ${PROJECT[]_HAS_IEEEFP_H} = true ; then
AC_DEFINE(PROJECT[]_HAS_IEEEFP_H,)
fi
])

dnl
dnl   check for ISO C99 stdint.h
dnl
AC_DEFUN(CF_CHECK_STDINT_H, [
AC_CHECK_HEADERS(stdint.h,
[PROJECT[]_HAS_STDINT_H=true],
[PROJECT[]_HAS_STDINT_H=false])
if test ${PROJECT[]_HAS_STDINT_H} = true ; then
AC_DEFINE(PROJECT[]_HAS_STDINT_H,)
fi
])

dnl
dnl   check whether values.h does really exist
dnl
AC_DEFUN(CF_CHECK_VALUES_H, [
AC_CHECK_HEADERS(values.h,
[PROJECT[]_HAS_VALUES_H=true],
[PROJECT[]_HAS_VALUES_H=false])
if test ${PROJECT[]_HAS_VALUES_H} = true ; then
AC_DEFINE(PROJECT[]_HAS_VALUES_H,)
fi
])

dnl
dnl   Check whether unistd.h does really exist.
dnl
AC_DEFUN(CF_CHECK_UNISTD_H, [
AC_CHECK_HEADERS(unistd.h,
[PROJECT[]_HAS_UNISTD_H=true],
[PROJECT[]_HAS_UNISTD_H=false])
if test ${PROJECT[]_HAS_UNISTD_H} = true ; then
AC_DEFINE(PROJECT[]_HAS_UNISTD_H,)
fi
])

dnl
dnl   Check whether limits.h does really exist.
dnl
AC_DEFUN(CF_CHECK_LIMITS_H, [
AC_CHECK_HEADERS(limits.h,
[PROJECT[]_HAS_LIMITS_H=true],
[PROJECT[]_HAS_LIMITS_H=false])
if test ${PROJECT[]_HAS_LIMITS_H} = true ; then
AC_DEFINE(PROJECT[]_HAS_LIMITS_H,)
fi
])

dnl
dnl   Check whether process.h does really exist.
dnl
AC_DEFUN(CF_CHECK_PROCESS_H, [
AC_CHECK_HEADERS(process.h,
[PROJECT[]_HAS_PROCESS_H=true],
[PROJECT[]_HAS_PROCESS_H=false])
if test ${PROJECT[]_HAS_PROCESS_H} = true ; then
AC_DEFINE(PROJECT[]_HAS_PROCESS_H,)
fi
])

dnl
dnl   Check whether sys/time.h does really exist.
dnl
AC_DEFUN(CF_CHECK_SYS_TIME_H, [
AC_CHECK_HEADERS(sys/time.h,
[PROJECT[]_HAS_SYS_TIME_H=true],
[PROJECT[]_HAS_SYS_TIME_H=false])
if test ${PROJECT[]_HAS_SYS_TIME_H} = true ; then
AC_DEFINE(PROJECT[]_HAS_SYS_TIME_H,)
fi
])

dnl
dnl   Check whether sys/stat.h does really exist.
dnl
AC_DEFUN(CF_CHECK_SYS_STAT_H, [
AC_CHECK_HEADERS(sys/stat.h,
[PROJECT[]_HAS_SYS_STAT_H=true],
[PROJECT[]_HAS_SYS_STAT_H=false])
if test ${PROJECT[]_HAS_SYS_STAT_H} = true ; then
AC_DEFINE(PROJECT[]_HAS_SYS_STAT_H,)
fi
])

dnl
dnl   Check whether sys/times.h does really exist.
dnl
AC_DEFUN(CF_CHECK_SYS_TIMES_H, [
AC_CHECK_HEADERS(sys/times.h,
[PROJECT[]_HAS_SYS_TIMES_H=true],
[PROJECT[]_HAS_SYS_TIMES_H=false])
if test ${PROJECT[]_HAS_SYS_TIMES_H} = true ; then
AC_DEFINE(PROJECT[]_HAS_SYS_TIMES_H,)
fi
])

dnl
dnl   Check whether sys/types.h does really exist.
dnl
AC_DEFUN(CF_CHECK_SYS_TYPES_H, [
AC_CHECK_HEADERS(sys/types.h,
[PROJECT[]_HAS_SYS_TYPES_H=true],
[PROJECT[]_HAS_SYS_TYPES_H=false])
if test ${PROJECT[]_HAS_SYS_TYPES_H} = true ; then
AC_DEFINE(PROJECT[]_HAS_SYS_TYPES_H,)
fi
])

dnl
dnl   Check whether sys/ioctl.h does really exist.
dnl
AC_DEFUN(CF_CHECK_SYS_IOCTL_H, [
AC_CHECK_HEADERS(sys/ioctl.h,
[PROJECT[]_HAS_SYS_IOCTL_H=true],
[PROJECT[]_HAS_SYS_IOCTL_H=false])
if test ${PROJECT[]_HAS_SYS_IOCTL_H} = true ; then
AC_DEFINE(PROJECT[]_HAS_SYS_IOCTL_H,)
fi
])

dnl
dnl   Check whether time.h does really exist.
dnl
AC_DEFUN(CF_CHECK_TIME_H, [
AC_CHECK_HEADERS(time.h,
[PROJECT[]_HAS_TIME_H=true],
[PROJECT[]_HAS_TIME_H=false])
if test ${PROJECT[]_HAS_TIME_H} = true ; then
AC_DEFINE(PROJECT[]_HAS_TIME_H,)
fi
])

dnl
dnl   Check whether sys/param.h does really exist.
dnl
AC_DEFUN(CF_CHECK_SYS_PARAM_H, [
AC_CHECK_HEADERS(sys/param.h,
[PROJECT[]_HAS_SYS_PARAM_H=true],
[PROJECT[]_HAS_SYS_PARAM_H=false])
if test ${PROJECT[]_HAS_SYS_PARAM_H} = true ; then
AC_DEFINE(PROJECT[]_HAS_SYS_PARAM_H,)
fi
])

dnl
dnl   Check whether dirent.h does really exist.
dnl
AC_DEFUN(CF_CHECK_DIRENT_H, [
AC_CHECK_HEADERS(dirent.h,
[PROJECT[]_HAS_DIRENT_H=true],
[PROJECT[]_HAS_DIRENT_H=false])
if test ${PROJECT[]_HAS_DIRENT_H} = true ; then
AC_DEFINE(PROJECT[]_HAS_DIRENT_H,)
fi
])

dnl
dnl   Check whether pwd.h does really exist.
dnl
AC_DEFUN(CF_CHECK_PWD_H, [
AC_CHECK_HEADERS(pwd.h,
[PROJECT[]_HAS_PWD_H=true],
[PROJECT[]_HAS_PWD_H=false])
if test ${PROJECT[]_HAS_PWD_H} = true ; then
AC_DEFINE(PROJECT[]_HAS_PWD_H,)
fi
])

dnl
dnl   Check whether direct.h does really exist.
dnl
AC_DEFUN(CF_CHECK_DIRECT_H, [
AC_CHECK_HEADERS(direct.h,
[PROJECT[]_HAS_DIRECT_H=true],
[PROJECT[]_HAS_DIRECT_H=false])
if test ${PROJECT[]_HAS_DIRECT_H} = true ; then
AC_DEFINE(PROJECT[]_HAS_DIRECT_H,)
fi
])

dnl
dnl   Check whether io.h does really exist.
dnl
AC_DEFUN(CF_CHECK_IO_H, [
AC_CHECK_HEADERS(io.h,
[PROJECT[]_HAS_IO_H=true],
[PROJECT[]_HAS_IO_H=false])
if test ${PROJECT[]_HAS_IO_H} = true ; then
AC_DEFINE(PROJECT[]_HAS_IO_H,)
fi
])

dnl
dnl   Check whether sys/socket.h does really exist.
dnl
AC_DEFUN(CF_CHECK_SYS_SOCKET_H, [
AC_CHECK_HEADERS(sys/socket.h,
[PROJECT[]_HAS_SYS_SOCKET_H=true],
[PROJECT[]_HAS_SYS_SOCKET_H=false])
if test ${PROJECT[]_HAS_SYS_SOCKET_H} = true ; then
AC_DEFINE(PROJECT[]_HAS_SYS_SOCKET_H,)
fi
])

dnl
dnl   Check whether netinet/in.h does really exist.
dnl
AC_DEFUN(CF_CHECK_NETINET_IN_H, [
AC_CHECK_HEADERS(netinet/in.h,
[PROJECT[]_HAS_NETINET_IN_H=true],
[PROJECT[]_HAS_NETINET_IN_H=false])
if test ${PROJECT[]_HAS_NETINET_IN_H} = true ; then
AC_DEFINE(PROJECT[]_HAS_NETINET_IN_H,)
fi
])

dnl
dnl   Check whether netdb.h does really exist.
dnl
AC_DEFUN(CF_CHECK_NETDB_H, [
AC_CHECK_HEADERS(netdb.h,
[PROJECT[]_HAS_NETDB_H=true],
[PROJECT[]_HAS_NETDB_H=false])
if test ${PROJECT[]_HAS_NETDB_H} = true ; then
AC_DEFINE(PROJECT[]_HAS_NETDB_H,)
fi
])

dnl
dnl   Check whether arpa/inet.h does really exist.
dnl
AC_DEFUN(CF_CHECK_ARPA_INET_H, [
AC_CHECK_HEADERS(arpa/inet.h,
[PROJECT[]_HAS_ARPA_INET_H=true],
[PROJECT[]_HAS_ARPA_INET_H=false])
if test ${PROJECT[]_HAS_ARPA_INET_H} = true ; then
AC_DEFINE(PROJECT[]_HAS_ARPA_INET_H,)
fi
])

dnl
dnl   Check whether sys/sysinfo.h does exist.
dnl
AC_DEFUN(CF_CHECK_SYS_SYSINFO_H, [
AC_CHECK_HEADERS(sys/sysinfo.h,
[PROJECT[]_HAS_SYS_SYSINFO_H=true],
[PROJECT[]_HAS_SYS_SYSINFO_H=false])
if test ${PROJECT[]_HAS_SYS_SYSINFO_H} = true ; then
AC_DEFINE(PROJECT[]_HAS_SYS_SYSINFO_H,)
fi
])

AC_DEFUN(CF_CHECK_SYSCONF, [
AC_CHECK_FUNCS(sysconf, HAS_SYSCONF=1)
if test "${HAS_SYSCONF}" = 1 ; then
AC_DEFINE(PROJECT[]_HAS_SYSCONF,)
fi
])

AC_DEFUN(CF_CHECK_KILL, [
AC_CHECK_FUNCS(kill, HAS_KILL=1)
if test "${HAS_KILL}" = 1 ; then
AC_DEFINE(PROJECT[]_HAS_KILL,)
fi
])

AC_DEFUN(CF_CHECK_HYPOT, [
AC_CHECK_FUNCS(kill, HAS_HYPOT=1)
if test "${HAS_HYPOT}" = 1 ; then
AC_DEFINE(PROJECT[]_HAS_HYPOT,)
fi
])

dnl
dnl   check whether vsnprintf is defined
dnl
AC_DEFUN(CF_CHECK_VSNPRINTF, [
AC_CHECK_FUNCS(vsnprintf, HAVE_VSNPRINTF=1)
if test "${HAVE_VSNPRINTF}" = 1 ; then
dnl
dnl   check whether vsnprintf works as expected
dnl   on Solaris 2.x it is broken in the 64bit version
dnl
AC_TRY_RUN(
[
#include <stdlib.h>
int main()
{
	char* buffer[50];
	vsnprintf(buffer, "%1s", "TEST");

	if (!strcmp(buffer, "T"))
	{
		return 0;
	}
	else
	{
		return 1;
	}
}
],
VSNPRINTF_OK=1,
DUMMY=0,
DUMMY=0
)

if test "${VSNPRINTF_OK}" = 1 ; then
AC_DEFINE(PROJECT[]_HAVE_VSNPRINTF)
fi
fi
])


dnl
dnl   check whether we need sysinfo or gethostname
dnl
AC_DEFUN(CF_CHECK_GETHOSTNAME, [
AC_CHECK_FUNCS(gethostname, HAVE_GETHOSTNAME=1)
if test "${HAVE_GETHOSTNAME}" = 1 ; then
AC_DEFINE(PROJECT[]_HAVE_GETHOSTNAME)
else
AC_CHECK_FUNCS(sysinfo, HAVE_SYSINFO=1)
if test "${HAVE_SYSINFO}" = 1  ; then
AC_DEFINE(PROJECT[]_HAVE_SYSINFO)
else
AC_MSG_RESULT()
AC_MSG_RESULT([Could not find gethostname or sysinfo methods!])
AC_MSG_RESULT([Please refer to config.log to identify the problem.])
CF_ERROR
fi
fi

dnl
dnl check for gethostname in the header
dnl
if test "${HAVE_GETHOSTNAME}" = 1 ; then
AC_MSG_CHECKING(for gethostname in unistd.h)
AC_TRY_COMPILE(
[
#include <unistd.h>
],
[
char name[1024];
gethostname(name, 1023);
],
HAVE_GETHOSTNAME_HEADER=1
)
if test "${HAVE_GETHOSTNAME_HEADER+set}" != set ; then
AC_MSG_RESULT(no)
AC_DEFINE(PROJECT[]_DEFINE_GETHOSTNAME)
else
AC_MSG_RESULT(yes)
fi
fi
])

AC_DEFUN(CF_CHECK_NETLIBS, [
dnl
dnl   first check if everythings already defined in libc
dnl
AC_CHECK_FUNCS(inet_addr, HAVE_INET_ADDR=1)
AC_CHECK_FUNCS(gethostbyname, HAVE_GETHOSTBYNAME=1)
if test "${HAVE_INET_ADDR+set}" = set ; then
AC_CHECK_FUNC(inet_aton, HAVE_INET_ATON=1)
fi

dnl   if gethostbyname was not defined in libc, try libxnet (Solaris only?)
if test "${HAVE_GETHOSTBYNAME+set}" != set -a "${USE_LIBXNET}" != false; then
AC_CHECK_LIB(xnet, gethostbyname)
unset ac_cv_func_gethostbyname
AC_CHECK_FUNCS(gethostbyname,HAVE_GETHOSTBYNAME=1)
fi
if test "${HAVE_INET_ADDR+set}" != set ; then
unset ac_cv_func_inet_addr
AC_CHECK_FUNCS(inet_addr,HAVE_INET_ADDR=1)
if test "${HAVE_INET_ADDR+set}" != set -a "${USE_LIBXNET}" != false; then
AC_CHECK_LIB(xnet, inet_addr)
unset ac_cv_func_inet_addr
AC_CHECK_FUNCS(inet_addr,HAVE_INET_ADDR=1)
fi
fi

if test "${HAVE_GETHOSTBYNAME+set}" != set ; then
AC_CHECK_LIB(nsl, gethostbyname)
unset ac_cv_func_gethostbyname
AC_CHECK_FUNCS(gethostbyname,HAVE_GETHOSTBYNAME=1)
fi

if test "${HAVE_INET_ADDR+set}" != set ; then
AC_CHECK_LIB(socket, inet_addr)
unset ac_cv_func_inet_addr
AC_CHECK_FUNCS(inet_addr,HAVE_INET_ADDR=1)
fi


dnl check again whether inet_aton exists (perhaps it was hidden in one
dnl of the other libraries..
if test "${HAVE_INET_ATON+set}" != set ; then
unset ac_cv_func_inet_aton
AC_CHECK_FUNC(inet_aton,HAVE_INET_ATON=1)
fi

if test "${HAVE_INET_ATON+set}" = set ; then
AC_DEFINE(HAVE_INET_ATON,)
fi

])

dnl
dnl   Check whether size arguments
dnl   in function calls like getsockname getpeername or accept
dnl   require a specialized typename or int
dnl   We simply compile a short example with all known types
dnl   and take one that didn't cause a warning (or an error)
dnl
AC_DEFUN(CF_CHECK_SOCKET_ARGS_AND_TYPES, [
AC_MSG_CHECKING(for socketlen type)
AC_TRY_COMPILE(
[
#include <sys/types.h>
#include <sys/socket.h>
#include <netinet/in.h>
],
[ sockaddr_in   addr;
socklen_t     len = 0;
getsockname(0, (struct sockaddr*)&addr, &len);
],
PROJECT[]_SOCKLEN_TYPE=socklen_t)
if test "${PROJECT[]_SOCKLEN_TYPE}" = "" ; then
AC_TRY_COMPILE(
[
#include <sys/types.h>
#include <sys/socket.h>
#include <netinet/in.h>
],
[ sockaddr_in   addr;
size_t    len = 0;
getsockname(0, (struct sockaddr*)&addr, &len);
],
PROJECT[]_SOCKLEN_TYPE=size_t)
fi
if test "${PROJECT[]_SOCKLEN_TYPE}" = "" ; then
AC_TRY_COMPILE(
[
#include <sys/types.h>
#include <sys/socket.h>
#include <netinet/in.h>
],
[ sockaddr_in   addr;
unsigned int  len = 0;
getsockname(0, (struct sockaddr*)&addr, &len);
],
PROJECT[]_SOCKLEN_TYPE="unsigned int")
fi
if test "${PROJECT[]_SOCKLEN_TYPE}" = "" ; then
AC_TRY_COMPILE(
[
#include <sys/types.h>
#include <sys/socket.h>
#include <netinet/in.h>
],
[ sockaddr_in   addr;
int           len = 0;
getsockname(0, (struct sockaddr*)&addr, &len);
],
PROJECT[]_SOCKLEN_TYPE="int")
fi
if test "${PROJECT[]_SOCKLEN_TYPE}" = "" ; then
AC_MSG_RESULT(FAILED)
AC_MSG_RESULT(-------------WARNING!---------------)
AC_MSG_RESULT(could not find a matching type for socket length argument)
AC_MSG_RESULT(in call to getsockname)
AC_MSG_RESULT(please check the setting for PROJECT[]_SOCKLEN_TYPE in config.mak)
AC_MSG_RESULT(and set it to the type needed for the third argument of getsockname)
AC_MSG_RESULT()
else
AC_MSG_RESULT($PROJECT[]_SOCKLEN_TYPE)
fi

AC_DEFINE_UNQUOTED(PROJECT[]_SOCKLEN_TYPE, ${PROJECT[]_SOCKLEN_TYPE})
])

dnl
dnl   check for the XDR functions: their interface and the libraries they're hidden in.
dnl
AC_DEFUN(CF_CHECK_XDR, [
if test "${PROJECT[]_NO_XDR}" = "true" ; then
AC_MSG_RESULT([No XDR headers available - building of XDR persistence support disabled])
AC_DEFINE(PROJECT[]_HAS_XDR, )
PROJECT[]_HAS_XDR=""
AC_SUBST(PROJECT[]_HAS_XDR)
else
AC_CHECK_HEADER(rpc/types.h, HAS_RPC_TYPES_H=true, HAS_RPC_TYPES_H=false)
if test "${HAS_RPC_TYPES_H}" = false ; then
AC_MSG_RESULT()
AC_MSG_RESULT([Cannot find RPC headers (rpc/types.h).])
AC_MSG_RESULT([If your system does not provide an RPC/XDR implementation (e.g., CYGWIN),])
AC_MSG_RESULT([please specify the option --without-xdr to avoid this error.])
CF_ERROR
fi

AC_CHECK_HEADER(rpc/xdr.h, HAS_XDR_H=true, HAS_XDR_H=false)
if test "${HAS_XDR_H}" = false ; then
AC_MSG_RESULT()
AC_MSG_RESULT([Cannot find XDR headers (rpc/xdr.h).])
CF_ERROR
fi

AC_MSG_CHECKING([arg types for xdrrec_create])
AC_TRY_COMPILE(
[
	#include <rpc/types.h>
	#include <rpc/xdr.h>
	extern "C" int dummy(void*, void*, unsigned int) {return 0;}
	void foo(){ 
		XDR xdrs;
		xdrrec_create(&xdrs, 0, 0, 0, dummy, dummy);
	}
],
[	
],
PROJECT[]_XDRREC_VOID_VOID_UINT=true,
PROJECT[]_XDRREC_VOID_VOID_UINT=false
)
if test "${PROJECT[]_XDRREC_VOID_VOID_UINT}" = true ; then
AC_MSG_RESULT([(void*, void*, unsigned int)])
else
AC_TRY_COMPILE(
	[
		 #include <rpc/types.h>
		 #include <rpc/xdr.h>
		 extern "C" int dummy(void*, char*, int) {return 0;}
		 void foo(){ 
			XDR xdrs;
			xdrrec_create(&xdrs, 0, 0, 0, dummy, dummy);
		 }
	 ],
	 [	
	 ],
	 PROJECT[]_XDRREC_VOID_CHAR_INT=true,
	 PROJECT[]_XDRREC_VOID_CHAR_INT=false
)
if test "${PROJECT[]_XDRREC_VOID_CHAR_INT}" = true ; then
	AC_MSG_RESULT([(void*, char*, int)])
else
	AC_TRY_COMPILE(
		[
			#include <rpc/types.h>
			#include <rpc/xdr.h>
			extern "C" int dummy(char*, char*, int) {return 0;}
			void foo(){ 
				XDR xdrs;
				xdrrec_create(&xdrs, 0, 0, 0, dummy, dummy);
			}
		],
		[	
		],
		PROJECT[]_XDRREC_CHAR_CHAR_INT=true,
		PROJECT[]_XDRREC_CHAR_CHAR_INT=false
	)
	if test "${PROJECT[]_XDRREC_CHAR_CHAR_INT}" = true ; then
		 AC_MSG_RESULT([(char*, char*, int)])
	else
		 AC_TRY_COMPILE(
			 [
				 #include <rpc/types.h>
				 #include <rpc/xdr.h>
				 extern "C" int dummy() {return 0;}
				 void foo(){
					XDR xdrs;
					xdrrec_create(&xdrs, 0, 0, 0, dummy, dummy);
				}
			 ],
			 [
			 ],
			 PROJECT[]_XDRREC_VOID=true,
			 PROJECT[]_XDRREC_VOID=false
		 )
		 if test "${PROJECT[]_XDRREC_VOID}" = true ; then
				AC_MSG_RESULT(())
		 else
	 AC_TRY_COMPILE(
		[
		 #include <rpc/types.h>
		 #include <rpc/xdr.h>
		extern "C" int dummy(void*, void*, int) {return 0;}
		void foo(){
		 XDR xdrs;
		 xdrrec_create(&xdrs, 0, 0, 0, dummy, dummy);
		}
	 ],
	 [
	 ],
	 PROJECT[]_XDRREC_VOID_VOID_INT=true,
	 PROJECT[]_XDRREC_VOID_VOID_INT=false
	 )
	if test "${PROJECT[]_XDRREC_VOID_VOID_INT}" = true ; then
		AC_MSG_RESULT((void*, void*, int))
	else
		AC_MSG_RESULT(not found!)
	fi
		fi 
fi
fi
fi

dnl
dnl  check whether there is a function to store 64bit
dnl  unsigned ints (xdr_u_hyper)
dnl
AC_MSG_CHECKING(for xdr_u_hyper function)
PROJECT[]_HAS_XDR_U_HYPER=false
AC_TRY_COMPILE(
[
#include <rpc/types.h>
#include <rpc/xdr.h>
],
[
xdr_u_hyper(0, 0);
],	
PROJECT[]_HAS_XDR_U_HYPER=true
)	

if test "${PROJECT[]_HAS_XDR_U_HYPER}" = "true" ; then
AC_MSG_RESULT(found)

AC_MSG_CHECKING([for 64-bit XDR type (for xdr_u_hyper)])
PROJECT[]_U_QUAD_TYPE=""
AC_TRY_COMPILE(
[
	#include <rpc/types.h>
	#include <rpc/xdr.h>
],
[	u_quad_t   q;
	XDR xdrs;
		xdr_u_hyper(&xdrs, &q);
],
PROJECT[]_U_QUAD_TYPE=u_quad_t
)	

if test "${PROJECT[]_U_QUAD_TYPE}" = "" ; then
AC_TRY_COMPILE(
	[
		#include <rpc/types.h>
		#include <rpc/xdr.h>
	],
	[	u_longlong_t   q;
		XDR xdrs;
		xdr_u_hyper(&xdrs, &q);
	],
	PROJECT[]_U_QUAD_TYPE=u_longlong_t
)	
fi

if test "${PROJECT[]_U_QUAD_TYPE}" = "" ; then
AC_TRY_COMPILE(
	[
		#include <rpc/types.h>
		#include <rpc/xdr.h>
	],
	[	unsigned long long int  q;
		XDR xdrs;
		xdr_u_hyper(&xdrs, &q);
	],
	PROJECT[]_U_QUAD_TYPE="unsigned long long int"
)	
fi

if test "${PROJECT[]_U_QUAD_TYPE}" = "" ; then
AC_TRY_COMPILE(
	[
		#include <rpc/types.h>
		#include <rpc/xdr.h>
	],
	[	__uint64_t  q;
		XDR xdrs;
		xdr_u_hyper(&xdrs, &q);
	],
	PROJECT[]_U_QUAD_TYPE=__uint64_t
)	
fi
if test "${PROJECT[]_U_QUAD_TYPE}" = "" ; then
AC_MSG_RESULT([Could not identify an appropriate type for xdr_u_hyper.])
CF_ERROR
fi

AC_MSG_RESULT(${PROJECT[]_U_QUAD_TYPE})
AC_DEFINE_UNQUOTED(PROJECT[]_XDR_UINT64_TYPE, ${PROJECT[]_U_QUAD_TYPE})
AC_DEFINE(PROJECT[]_HAS_XDR_U_HYPER)

else

dnl
dnl	we do not have xdr_u_hyper, so PROJECT[] has to use two 
dnl	calls to xdr_u_int instead. 
dnl	However, we have to identify whether the system supports
dnl	64bit unsigned types at all
dnl	
AC_TRY_COMPILE(
[
],
[	
	unsigned long long int  q = 1234567890123456789LL;
],
PROJECT[]_U_QUAD_TYPE="unsigned long long int"
)	
if test "${PROJECT[]_U_QUAD_TYPE}" = "" ; then
AC_MSG_RESULT([Could not identify an 64 bit unsigned type (long long).])
CF_ERROR
fi

AC_DEFINE_UNQUOTED(PROJECT[]_XDR_UINT64_TYPE, ${PROJECT[]_U_QUAD_TYPE})
AC_MSG_RESULT(unsigned long long int)

fi


dnl
dnl Define appropriate symbols in config.h.
dnl The symbols are used in CONCEPT/XDRPersistenceManager.C only.
dnl
if test "${PROJECT[]_XDRREC_VOID_CHAR_INT}" = true ; then
AC_DEFINE(PROJECT[]_XDRREC_CREATE_VOID_CHAR_INT)
fi
if test "${PROJECT[]_XDRREC_CHAR_CHAR_INT}" = true ; then
AC_DEFINE(PROJECT[]_XDRREC_CREATE_CHAR_CHAR_INT)
fi
if test "${PROJECT[]_XDRREC_VOID}" = true ; then
AC_DEFINE(PROJECT[]_XDRREC_CREATE_VOID)
fi
if test "${PROJECT[]_XDRREC_VOID_VOID_UINT}" = true ; then
AC_DEFINE(PROJECT[]_XDRREC_CREATE_VOID_VOID_UINT)
fi
if test "${PROJECT[]_XDRREC_VOID_VOID_INT}" = true ; then
AC_DEFINE(PROJECT[]_XDRREC_CREATE_VOID_VOID_INT)
fi

dnl
dnl		Try to guess the library required for the XDR stuff.
dnl   It is often in libc, but Solaris hides it in libnsl.
dnl
AC_MSG_CHECKING(for XDR symbols in libc)
AC_TRY_LINK([
#include <rpc/types.h>
#include <rpc/xdr.h>
],
[
XDR xdrs;
int i;
xdr_int(&xdrs, &i);
],
XDR_IN_LIBC=true
)
if test "${XDR_IN_LIBC}" != true ; then
AC_MSG_RESULT(not found!)
AC_MSG_CHECKING(for XDR symbols in libnsl)

SAVE_LIBS=${LIBS}
LIBS="-lnsl ${LIBS}"
AC_TRY_LINK([
	#include <rpc/types.h>
	#include <rpc/xdr.h>
],
[
	XDR xdrs;
	int i;
	xdr_int(&xdrs, &i);
],
XDR_IN_LIBNSL=true
)
if test "${XDR_IN_LIBNSL}" = true ; then
AC_MSG_RESULT(yes)
else
AC_MSG_RESULT(no)
AC_MSG_RESULT()
AC_MSG_RESULT(Did not find XDR symbols in libc or libnsl.)
CF_ERROR
fi
else
dnl
dnl	XDR symbols are in libc.
dnl
AC_MSG_RESULT(yes)
fi
AC_DEFINE(PROJECT[]_HAS_XDR, true)
PROJECT[]_HAS_XDR=true
AC_SUBST(PROJECT[]_HAS_XDR)
fi
])

dnl
dnl FFTW -- Fastest Fourier Transform in the West
dnl
AC_DEFUN(CF_CHECK_FFTW_SUPPORT, [
if test "${FFTW_SUPPORT}" = true ; then
	FFTW_SUPPORT=true

	AC_MSG_CHECKING(for FFTW headers)
	if test "${FFTW_INCL}" != "" ; then
		CF_FIND_HEADER(FFTW_INCL_PATH, fftw3.h, ${FFTW_INCL})
	else
		CF_FIND_HEADER(FFTW_INCL_PATH, fftw3.h,)
	fi

	if test "${FFTW_INCL_PATH}" = "" ; then
		AC_MSG_RESULT((not found!))
		AC_MSG_RESULT()
		AC_MSG_RESULT(No FFTW header files found! Please specify the path to the FFTW header fftw3.h)
		AC_MSG_RESULT(by passing the option --with-fftw-incl=DIR to configure.)
		AC_MSG_RESULT(The FFTW package can be found under the following URL:)
		AC_MSG_RESULT(  http://www.fftw.org)
		CF_ERROR
	else
		AC_MSG_RESULT((${FFTW_INCL_PATH}))
	fi

	if test "${FFTW_DISABLE_FFTW_FLOAT}" = "false" ; then
		AC_MSG_CHECKING(for FFTW library with float support)
		CF_FIND_LIB(FFTW_LIB_F,libfftw3f, ${FFTW_LIBPATH})

		if test "${FFTW_LIB_F}" = "" ; then
			AC_MSG_RESULT((not found!))
			FFTW_DISABLE_FFTW_FLOAT=true
		else
			AC_MSG_RESULT((${FFTW_LIB_F}))
			FFTW_DISABLE_FFTW_FLOAT=false
		fi
	fi

	if test "${FFTW_DISABLE_FFTW_DOUBLE}" = "false" ; then
		AC_MSG_CHECKING(for FFTW library with double support)
		CF_FIND_LIB(FFTW_LIB_D,libfftw3, ${FFTW_LIBPATH})

		if test "${FFTW_LIB_D}" = "" ; then
			AC_MSG_RESULT((not found!))
			FFTW_DISABLE_FFTW_DOUBLE=true
		else
			AC_MSG_RESULT((${FFTW_LIB_D}))
			FFTW_DISABLE_FFTW_DOUBLE=false
		fi
	fi

	if test "${FFTW_DISABLE_FFTW_LONGDBL}" = false ; then
		AC_MSG_CHECKING(for FFTW library with long double support)
		CF_FIND_LIB(FFTW_LIB_L,libfftw3l, ${FFTW_LIBPATH})

		if test "${FFTW_LIB_L}" = "" ; then
			AC_MSG_RESULT((not found!))
			FFTW_DISABLE_FFTW_LONGDBL=true
		else
			AC_MSG_RESULT((${FFTW_LIB_L}))
			FFTW_DISABLE_FFTW_LONGDBL=false
		fi
	fi

	if test "${FFTW_LIB_F}" = "" -a "${FFTW_LIB_D}" = "" -a "${FFTW_LIB_L}" = "" ; then
			AC_MSG_RESULT(Please install it in a standard location or specify the path)
			AC_MSG_RESULT(with --with-fftw-lib=DIR on the command line.)
			CF_ERROR
	fi

	dnl prevent the use of -L/usr/lib - this may lead to problems with different
	dnl binary formats (e.g. SGI O32/N32 format)
	if test "${FFTW_INCL_PATH}" != /usr/include -a "${FFTW_INCL_PATH}" != "" ; then
		PROJECT[]_INCLUDES="${PROJECT[]_INCLUDES} -I${FFTW_INCL_PATH}"
	fi

	AC_DEFINE(PROJECT[]_HAS_FFTW, true)
	AC_DEFINE(PROJECT[]_HAS_FFTW_H, true)
	PROJECT[]_HAS_FFTW=true
	AC_SUBST(PROJECT[]_HAS_FFTW)
	AC_SUBST(PROJECT[]_HAS_FFTW_H)

	PROJECT[]_HAS_FFTW_FLOAT=""
	if test "${FFTW_DISABLE_FFTW_FLOAT}" = "false" ; then
		AC_MSG_CHECKING(linking against libfftw3f)
		SAVE_LIBS=${LIBS}
		SAVE_LDFLAGS=${LDFLAGS}
		LIBS="${FFTW_LIB_F}/libfftw3f.a ${LIBS}"
		LDFLAGS="$LDFLAGS -I${FFTW_INCL_PATH}"
		FFTW_LINKING_OK=0
		AC_TRY_LINK([
		#include <fftw3.h>
		],
		[
		 fftwf_plan f = fftwf_plan_dft_1d(1,0,0,1,FFTW_FORWARD);
		], FFTW_LINKING_OK=1)
		LIBS=${SAVE_LIBS}
		LDFLAGS=${SAVE_LDFLAGS}
		if test "${FFTW_LINKING_OK}" != "1" ; then
	AC_MSG_RESULT(no)
	AC_MSG_RESULT()
	AC_MSG_RESULT([Cannot link against libfftw3f. Please check config.log and])
	AC_MSG_RESULT([specify appropriate options to configure (e.g. --with-fftw-lib/incl).])
	CF_ERROR
	else
	AC_MSG_RESULT(yes)
	PROJECT[]_HAS_FFTW_FLOAT=true
	fi
	fi

	PROJECT[]_HAS_FFTW_DOUBLE=""
	if test "${FFTW_DISABLE_FFTW_DOUBLE}" = "false" ; then
	AC_MSG_CHECKING(linking against libfftw3)
	SAVE_LIBS=${LIBS}
	SAVE_LDFLAGS=${LDFLAGS}
	LIBS="${FFTW_LIB_D}/libfftw3.a ${LIBS}"
	LDFLAGS="$LDFLAGS -I${FFTW_INCL_PATH}"
	FFTW_LINKING_OK=0
	AC_TRY_LINK([
							#include <fftw3.h>
						],
						[
							 fftw_plan f = fftw_plan_dft_1d(1,0,0,1,FFTW_FORWARD);
						], FFTW_LINKING_OK=1)
	LIBS=${SAVE_LIBS}
	LDFLAGS=${SAVE_LDFLAGS}
	if test "${FFTW_LINKING_OK}" != "1" ; then
	AC_MSG_RESULT(no)
	AC_MSG_RESULT()
	AC_MSG_RESULT([Cannot link against libfftw3. Please check config.log and])
	AC_MSG_RESULT([specify appropriate options to configure (e.g. --with-fftw-lib/incl).])
	CF_ERROR
	else
	AC_MSG_RESULT(yes)
	PROJECT[]_HAS_FFTW_DOUBLE=true
	fi
	fi

	PROJECT[]_HAS_FFTW_LONG_DOUBLE=""
	if test "${FFTW_DISABLE_FFTW_LONGDBL}" = "false" ; then
	AC_MSG_CHECKING(linking against libfftw3l)
	SAVE_LIBS=${LIBS}
	SAVE_LDFLAGS=${LDFLAGS}
	LIBS="${FFTW_LIB_D}/libfftw3l.a ${LIBS}"
	LDFLAGS="$LDFLAGS -I${FFTW_INCL_PATH}"
	FFTW_LINKING_OK=0
	AC_TRY_LINK([
							#include <fftw3.h>
						],
						[
							 fftwl_plan f = fftwl_plan_dft_1d(1,0,0,1,FFTW_FORWARD);
						], FFTW_LINKING_OK=1)
	LIBS=${SAVE_LIBS}
	LDFLAGS=${SAVE_LDFLAGS}
	if test "${FFTW_LINKING_OK}" != "1" ; then
	AC_MSG_RESULT(no)
	AC_MSG_RESULT()
	AC_MSG_RESULT([Cannot link against libfftw3l. Please check config.log and])
	AC_MSG_RESULT([specify appropriate options to configure (e.g. --with-fftw-lib/incl).])
	CF_ERROR
	else
	AC_MSG_RESULT(yes)
	PROJECT[]_HAS_FFTW_LONG_DOUBLE=true
	fi
	fi
fi
AC_DEFINE_UNQUOTED(PROJECT[]_COMPLEX_TYPE, ${PROJECT[]_COMPLEX_TYPE})
if test "${PROJECT[]_HAS_FFTW_FLOAT}" != "" ; then
FFTW_LIBS="${FFTW_LIB_F}/libfftw3f.a}"
AC_DEFINE(PROJECT[]_HAS_FFTW_FLOAT,)
fi
if test "${PROJECT[]_HAS_FFTW_DOUBLE}" != "" ; then
FFTW_LIBS="${FFTW_LIB_D}/libfftw3.a"
AC_DEFINE(PROJECT[]_HAS_FFTW_DOUBLE,)
fi
if test "${PROJECT[]_HAS_FFTW_LONG_DOUBLE}" != "" ; then
FFTW_LIBS="${FFTW_LIB_L}/libfftw3l.a"
AC_DEFINE(PROJECT[]_HAS_FFTW_LONG_DOUBLE,)
fi
if test "${PROJECT[]_HAS_FFTW_H}" != "" ; then
LDFLAGS="$LDFLAGS -I${FFTW_INCL_PATH}"
fi

AC_SUBST(PROJECT[]_HAS_FFTW_FLOAT,)
AC_SUBST(PROJECT[]_HAS_FFTW_DOUBLE,)
AC_SUBST(PROJECT[]_HAS_FFTW_LONG_DOUBLE,)
AC_SUBST(FFTW_LIBS)
])

dnl
dnl		VIEW support
dnl
AC_DEFUN(CF_VIEW, [
dnl
dnl    search for X-libs and includes and PROJECT[]View (OpenGL/MESA) stuff
dnl 
if test "${USE_VIEW}" = true ; then
AC_PATH_X
X11_INCPATH=${x_includes}
X11_LIBPATH=${x_libraries}

if test "${no_x}" = "yes" ; then
USE_VIEW=false
fi

if test "${USE_VIEW}" = true ; then
if test "${X11_LIBPATH}" = "/usr/lib" -o "${X11_LIBPATH}" = "" ; then
X11_LIBPATH=""
X11_LIBPATHOPT=""
else
X11_LIBPATHOPT="-L${X11_LIBPATH}"
fi
fi

if test "${USE_VIEW}" = true ; then
dnl
dnl		Fix up the OpenGL stuff for MacOS X -- here we need to use OpenGL and AGL frameworks
dnl
if test "${OS}" = "Darwin" ; then
	VIEW_PLATFORM="OpenGL-Darwin"
	OPENGL_LIBOPTS="-framework Carbon -framework OpenGL -framework AGL"
	X11_LIBPATHOPT=""
	X11_INCPATH=""
	OPENGL_INCPATH="-I/System/Library/Frameworks/OpenGL.framework/Versions/A/Headers -I/System/Library/Frameworks/AGL.framework/Versions/A/Headers"
	VIEW_INCLUDES="${VIEW_INCLUDES} ${OPENGL_INCPATH}"
fi

if test "${VIEW_PLATFORM}" = Mesa ; then
AC_MSG_CHECKING(for Mesa includes)
CF_FIND_HEADER(MESA_INCLUDES,GL/gl.h, ${OPENGL_INCPATH} ${X11_INCPATH})
if test "${MESA_INCLUDES}" = "" ; then
	AC_MSG_RESULT((not found!))
	AC_MSG_RESULT()
	AC_MSG_RESULT(No Mesa headers found! Please specify the path to the directory)
	AC_MSG_RESULT(containing the Mesa headers using --with-opengl-incl=DIR.)
	AC_MSG_RESULT(Mesa can be obtained from www.mesa3d.org.)
	CF_ERROR
else
	AC_MSG_RESULT(${MESA_INCLUDES})
fi

if test "${USE_VIEW}" = true ; then
	AC_MSG_CHECKING(for Mesa library)
	CF_FIND_LIB(MESA_LIBS,libMesaGL, ${OPENGL_LIBPATH} ${X11_LIBPATH})
	if test "${MESA_LIBS}" = "" ; then
		CF_FIND_LIB(MESA_LIBS,libGL, ${OPENGL_LIBPATH} ${X11_LIBPATH})
	fi
	if test "${MESA_LIBS}" = "" ; then
		AC_MSG_RESULT((not found!))
		AC_MSG_RESULT()
		AC_MSG_RESULT(No Mesa library libMesaGL or libGL found! Please specify the path)
		AC_MSG_RESULT(to the directory containing the library using the --with-opengl-libs=DIR.)
		AC_MSG_RESULT(Mesa can be obtained from www.mesa3d.org.)
		AC_MSG_RESULT(Aborted.)
	else
		AC_MSG_RESULT((${MESA_LIBS}))
	fi
fi

dnl prevent the use of -L/usr/lib - this may lead to problems with different
dnl binary formats (e.g. SGI O32/N32 format)
if test "${MESA_INCLUDES}" != /usr/include -a "${MESA_INCLUDES}" != "" ; then
	VIEW_INCLUDES="${VIEW_INCLUDES} -I${MESA_INCLUDES}"
fi
fi

if test ${VIEW_PLATFORM} = OpenGL ; then
AC_MSG_CHECKING(for OpenGL includes)
CF_FIND_HEADER(OPENGL_INCPATH,GL/gl.h)
if test "${OPENGL_INCPATH}" = "" ; then
	AC_MSG_RESULT((not found!))
	AC_MSG_RESULT()
	AC_MSG_RESULT(no OpenGL headers found! Please use the option --with-opengl-incl=DIR)
	AC_MSG_RESULT(of configure to specify the correct path to these headers.)
	CF_ERROR
else
	AC_MSG_RESULT((${OPENGL_INCPATH}))
fi

AC_MSG_CHECKING(for OpenGL library)
CF_FIND_LIB(OPENGL_LIBPATH,libGL)
if test "${OPENGL_LIBPATH}" = "" ; then
	AC_MSG_RESULT((not found!))
	AC_MSG_RESULT()
	AC_MSG_RESULT(no OpenGL lib found! Please use the option --with-opengl-libs=DIR)
	AC_MSG_RESULT(of configure to specify the correct path to these libraries.)
	CF_ERROR
else
	AC_MSG_RESULT((${OPENGL_LIBPATH}))
fi

if test "${OPENGL_INCPATH}" != /usr/include && test "${OPENGL_INCPATH}" != "" ; then
	VIEW_INCLUDES="${VIEW_INCLUDES} -I${OPENGL_INCPATH}"
fi
fi

if test "${USE_VIEW}" = true ; then
AC_MSG_CHECKING(for QT headers)
if test "${QTDIR}" != "" ; then
	CF_FIND_HEADER(QT_INCPATH,qgl.h,${QTDIR}/include ${PROJECT[]_PATH}/contrib/qt/include)
else
	CF_FIND_HEADER(QT_INCPATH,qgl.h,${PROJECT[]_PATH}/contrib/qt/include)
fi

if test "${QT_INCPATH}" = "" ; then
	AC_MSG_RESULT((not found!))
	AC_MSG_RESULT()
	AC_MSG_RESULT(No QT header files found! Please specify the path to the QT headers)
	AC_MSG_RESULT(by passing the option --with-qt-incl=DIR to configure.)
	AC_MSG_RESULT(You may also set the environment variable QTDIR to the correct)
	AC_MSG_RESULT(path - configure will recognize this, too.)
	AC_MSG_RESULT(The QT package can be found under the following URL:)
	AC_MSG_RESULT(  http://www.troll.no/qt)
	CF_ERROR
else
	AC_MSG_RESULT((${QT_INCPATH}))	
fi

AC_MSG_CHECKING(for libqt${QT_MT_SUFFIX})
if test "${QTDIR}" != "" ; then
	if test "${QT_MT_SUFFIX}" = "-mt" -o "${QT_MT_SUFFIX+set}" != set ; then
						if test -a "${QTDIR}/lib/libqt-mt.so" ; then
							QT_LIBPATH="${QTDIR}/lib"
						elif test -a "${QTDIR}/lib/libqt.so" ; then
							QT_LIBPATH="${QTDIR}/lib"
							QT_MT_SUFFIX=""
						fi
						if test "${QT_LIBPATH}" = "" ; then
							CF_FIND_LIB(QT_LIBPATH, libqt${QT_MT_SUFFIX}, ${QTDIR}/lib ${QTDIR}/lib/${BINFMT} ${PROJECT[]_PATH}/contrib/qt/include)
						fi
					else	
						if test -a "${QTDIR}/lib/libqt.so" ; then
							QT_LIBPATH="${QTDIR}/lib"
						elif test -a "${QTDIR}/lib/libqt-mt.so" ; then
							QT_LIBPATH="${QTDIR}/lib"
							QT_MT_SUFFIX="-mt"
						fi
						if test "${QT_LIBPATH}" = "" ; then
							CF_FIND_LIB(QT_LIBPATH, libqt${QT_MT_SUFFIX}, ${QTDIR}/lib ${QTDIR}/lib/${BINFMT} ${PROJECT[]_PATH}/contrib/qt/include)
						fi
					fi
				else
					CF_FIND_LIB(QT_LIBPATH, libqt${QT_MT_SUFFIX}, ${PROJECT[]_PATH}/contrib/qt/lib ${PROJECT[]_PATH}/contrib/qt/lib/${BINFMT})
				fi

				if test "${QT_LIBPATH}" = "" ; then
					AC_MSG_RESULT((not found!))
					AC_MSG_RESULT()
					AC_MSG_RESULT([The QT library could not be found. Please specify the path to libqt])
					AC_MSG_RESULT([by passing the option --with-qt-libs=DIR to configure.])
					AC_MSG_RESULT([You may also set the environment variable QTDIR to the correct])
					AC_MSG_RESULT([path - configure will recognize this, too.])
					AC_MSG_RESULT([If the QT library was built with thread support enabled (libqt-mt])
					AC_MSG_RESULT([instead of libqt), please specify the option --with-threadsafe-qt.])
					AC_MSG_RESULT([The QT package can be found under the following URL:])
					AC_MSG_RESULT(  http://www.troll.no/qt)
					CF_ERROR
				else
					AC_MSG_RESULT((${QT_LIBPATH}))	
				fi

				
				dnl
				dnl extract the QT version number and version number string from include/qglobal.h
				dnl
				QT_VERSION=`${GREP} "#define QT_VERSION[^_]" ${QT_INCPATH}/qglobal.h | ${TR} '\011' ' ' | ${TR} -s ' ' | ${CUT} -d\  -f3`
				QT_VERSION_STR=`${GREP} "#define QT_VERSION_STR" ${QT_INCPATH}/qglobal.h | ${TR} '\011' ' ' | ${TR} -s ' ' | ${CUT} -d\  -f3 | ${TR} -d '\042'`
				AC_MSG_CHECKING(for QT version number in qglobal.h)
				if test "${QT_VERSION}" = "" ; then
					AC_MSG_RESULT([<unknown>])
					AC_MSG_RESULT()
					AC_MSG_RESULT([  Could not determine version number of QT library -- please])
					AC_MSG_RESULT([  check config.log for details.])
					AC_MSG_RESULT([  You might have a problem with your (DY)LD_LIBRARY_PATH.])
					AC_MSG_RESULT([  Please check the settings of QTDIR as well or specify])
					AC_MSG_RESULT([  the path to the library/headers with])
					AC_MSG_RESULT([    --with-qt-libs=<DIR> / --with-qt-incl=<DIR>])
					CF_ERROR
				else
					AC_MSG_RESULT([${QT_VERSION} (${QT_VERSION_STR})])
					AC_DEFINE_UNQUOTED(PROJECT[]_QT_VERSION, ${QT_VERSION})
					AC_DEFINE_UNQUOTED(PROJECT[]_QT_VERSION_STR, ${QT_VERSION_STR})
					if test "${QT_MT_SUFFIX}" = "-mt" ; then
						AC_DEFINE(PROJECT[]_QT_HAS_THREADS,)
					fi
				fi			
		
				dnl
				dnl  We do require QT 3.x by now. 2.x won't do...
				dnl
				if test `echo ${QT_VERSION} | ${CUT} -c1-2` != "0x" ; then
					if test "${QT_VERSION}" -lt 300 ; then
						AC_MSG_RESULT()
						AC_MSG_RESULT([QT version 3.0 or above is required for PROJECT[]. Please update])
						AC_MSG_RESULT([to a more current version or specify the path to a more])
						AC_MSG_RESULT([recent version of libqt by passing the option --with-qt-libs=DIR])
						AC_MSG_RESULT([to configure.])
						AC_MSG_RESULT([You may also set the environment variable QTDIR to the correct])
						AC_MSG_RESULT([path - configure will recognize this, too.])
						AC_MSG_RESULT()
						AC_MSG_RESULT([The complete QT package can be found under the following URL:])
						AC_MSG_RESULT([  http://www.troll.no/qt])
						CF_ERROR
					fi
				fi

				dnl
				dnl	Add the QT include path to the VIEW includes
				dnl
				if test "${QT_INCPATH}" != /usr/include && test "${QT_INCPATH}" != "" ; then
					VIEW_INCLUDES="${VIEW_INCLUDES} -I${QT_INCPATH}"
				fi	
			fi
		fi
	fi


	dnl
	dnl   verify libraries needed for VIEW
	dnl   (X, QT, Mesa/OpenGL)
	dnl

	if test "${USE_VIEW}" = true ; then		
		dnl  
		dnl
		dnl  identify the X11 libraries needed to link agains
		dnl
		dnl
		
		AC_MSG_CHECKING(linking against X11 libraries)
		dnl 
		dnl   if the user specified X libraries, try these first
		dnl
		if test "${X11_LIBS}" != "" ; then
			SAVE_LIBS=${LIBS}
			LIBS="${X11_LIBPATHOPT} ${X11_LIBS} ${LIBS}"
			AC_TRY_LINK([],[],X_LINKING_OK=1)
			LIBS=${SAVE_LIBS}
		fi

		dnl
		dnl  Special treatment for MacOS X -- we just ignore everything.
		dnl		The OpenGL and AGL frameworks will take care of it...
		dnl
		if test "${OS}" = "Darwin" ; then
			X11_LIBS=""
			X11_LIBPATHOPTS=""
			X_LINKING_OK="true"
		fi

		dnl 		
		dnl  now try the default guess: Xmu, Xext, Xt, and X11 
		dnl
		if test "${X_LINKING_OK+set}" != set ; then
			X11_LIBS="-lXmu -lXext -lXt -lX11 -lm"
			SAVE_LIBS=${LIBS}
			LIBS="${X11_LIBPATHOPT} ${X11_LIBS} ${LIBS}"
			AC_TRY_LINK([],[],X_LINKING_OK=1)
			LIBS=${SAVE_LIBS}
		fi
		
		dnl 		
		dnl  second guess: add SM and ICE
		dnl
		if test "${X_LINKING_OK+set}" != set ; then
			X11_LIBS="-lXmu -lXext -lXt -lX11 -lSM -lICE -lm"
			SAVE_LIBS=${LIBS}
			LIBS="${X11_LIBPATHOPT} ${X11_LIBS} ${LIBS}"
			AC_TRY_LINK([],[],X_LINKING_OK=1)
			LIBS=${SAVE_LIBS}
		fi
		
		dnl 		
		dnl  now try the default guess: Xmu, Xext, Xt, and X11 
		dnl
		if test "${X_LINKING_OK+set}" != set ; then
			X11_LIBS="-lXmu -lXt -lX11 -lm"
			SAVE_LIBS=${LIBS}
			LIBS="${X11_LIBPATHOPT} ${X11_LIBS} ${LIBS}"
			AC_TRY_LINK([],[],X_LINKING_OK=1)
			LIBS=${SAVE_LIBS}
		fi
		
		dnl 		
		dnl  second guess: add SM and ICE
		dnl
		if test "${X_LINKING_OK+set}" != set ; then
			X11_LIBS="-lXmu -lXt -lX11 -lSM -lICE -lm"
			SAVE_LIBS=${LIBS}
			LIBS="${X11_LIBPATHOPT} ${X11_LIBS} ${LIBS}"
			AC_TRY_LINK([],[],X_LINKING_OK=1)
			LIBS=${SAVE_LIBS}
		fi
		
		dnl 
		dnl  if we could not link - complain about it!
		dnl
		if test "${X_LINKING_OK+set}" = set ; then
			AC_MSG_RESULT(yes)	
		else
			AC_MSG_RESULT(no)
			AC_MSG_RESULT()
			AC_MSG_RESULT(Don't know how to link with X11 libraries.)
			AC_MSG_RESULT(Please specify the correct libraries (e.g. -lXmu -lXt -lX11) in the)
			AC_MSG_RESULT(environment variable X11_LIBS)
			AC_MSG_RESULT(If you are running Solaris 2.x you might also try the option --without-libxnet)
			AC_MSG_RESULT(if your X libraries were linked against libsocket and libnsl instead of libxnet.)
			AC_MSG_RESULT(Built of visualization component VIEW disabled.)
			AC_MSG_RESULT()
			USE_VIEW=false
		fi

		dnl		
		dnl  define some variables: X11_LIBOPTS and VIEW_LIBS
		dnl
		X11_LIBOPTS="${X11_LIBPATHOPT} ${X11_LIBS}"
	fi

	if test "${USE_VIEW}" = true ; then
		if test "${VIEW_PLATFORM}" = OpenGL ; then
			if test "${OPENGL_LIBPATH}" != "/usr/lib" -a "${OPENGL_LIBPATH}" != "" ; then
				OPENGL_LIBOPTS="-L${OPENGL_LIBPATH} -lGLU -lGL"
			else
				OPENGL_LIBPATH=""
				OPENGL_LIBOPTS="-lGLU -lGL"
			fi

			dnl make sure we have OpenGL libs and no Mesa libs!
			dnl
			SAVE_LIBS=${LIBS}
			SAVE_LDFLAGS=${LDFLAGS}
			LIBS="${LIBS} ${X11_LIBOPTS}"
			if test "${OPENGL_LIBPATH}" != "" ; then
				LDFLAGS="${LDFLAGS} -L${OPENGL_LIBPATH}"
			fi
			AC_CHECK_LIB(GL, XMesaGarbageCollect, VIEW_PLATFORM=Mesa)
			LIBS=${SAVE_LIBS}
			LDFLAGS=${SAVE_LDFLAGS}
			if test "${VIEW_PLATFORM}" != Mesa ; then
				AC_MSG_CHECKING(linking against OpenGL libraries)
				SAVE_LIBS=${LIBS}
				LIBS="${OPENGL_LIBOPTS} ${LIBS}"
				AC_TRY_LINK([],[],OPENGL_LINKING_OK=1)
				LIBS=${SAVE_LIBS}
				if test "${OPENGL_LINKING_OK+set}" != set ; then
					AC_MSG_RESULT(no)
					AC_MSG_RESULT()
					AC_MSG_RESULT(Cannot link against libGL/GLU - disabling visualization support!)
					AC_MSG_RESULT(Please specify the path to OpenGL libraries using --with-opengl-libs=DIR)
					CF_ERROR
				else
					AC_MSG_RESULT(yes)
				fi
			fi
		fi
	fi

	if test "${USE_VIEW}" = true ; then
		if test "${VIEW_PLATFORM}" = Mesa ; then
			dnl
			dnl  strip default path
			dnl
		
			if test "${MESA_LIBS}" = "" ; then 
				MESA_LIBS=${OPENGL_LIBPATH}
			fi
			if test "${MESA_LIBS}" != "/usr/lib" -a "${MESA_LIBS}" != "" ; then
				OPENGL_LIBPATH="${MESA_LIBS}"
				OPENGL_LIBPATHOPT="-L${MESA_LIBS}"			
			else
				OPENGL_LIBPATH=""
				OPENGL_LIBPATHOPT=""
			fi
			
			dnl
			dnl  out first guess for the names of the Mesa libraries
			dnl
			OPENGL_LIBS="-lGLU -lGL"

			dnl
			dnl  try to link against mesa libraries
			dnl
			AC_MSG_CHECKING(linking against Mesa libs)
			SAVE_LIBS=${LIBS}
			LIBS="${OPENGL_LIBPATHOPT} ${OPENGL_LIBS} ${X11_LIBOPTS} ${LIBS} "
			AC_TRY_LINK([],[], HAVE_MESALIBS=1)
			LIBS=${SAVE_LIBS}

			dnl
			dnl  could not link against libGLU/libGL,
			dnl  so try libMesaGLU/libMesaGL
			dnl
			if test "${HAVE_MESALIBS+set}" != set ; then
				OPENGL_LIBS="-lMesaGLU -lMesaGL"
				SAVE_LIBS=${LIBS}
				LIBS="${OPENGL_LIBPATHOPT} ${OPENGL_LIBS} ${X11_LIBOPTS} ${LIBS} "
				AC_TRY_LINK([],[], HAVE_MESALIBS=1)
				LIBS=${SAVE_LIBS}
			fi

			if test "${HAVE_MESALIBS+set}" != set ; then
				AC_MSG_RESULT(no)
				AC_MSG_RESULT()
				AC_MSG_RESULT(Cannot link against libMesaGL/GLU - disabling visualization support!)
				AC_MSG_RESULT(Please specify the path to libMesaGL using --with-opengl-libs=DIR)
				CF_ERROR
			else
				AC_MSG_RESULT(yes)
				OPENGL_LIBOPTS="${OPENGL_LIBPATHOPT} ${OPENGL_LIBS}"
			fi
		fi
	fi

	if test "${USE_VIEW}" = true ; then
		if test "${QT_LIBPATH}" != "/usr/lib" ; then
			QTQGL_LIBOPTS="-L${QT_LIBPATH} -lqgl -lqt${QT_MT_SUFFIX}"
			QT_LIBOPTS="-L${QT_LIBPATH} -lqt${QT_MT_SUFFIX}"
		else 
			QT_LIBPATH=""
			QTQGL_LIBOPTS="-lqgl -lqt${QT_MT_SUFFIX}"
			QT_LIBOPTS="-lqt${QT_MT_SUFFIX}"
		fi
	fi

	if test "${USE_VIEW}" = true ; then
		AC_MSG_CHECKING(linking against QT libraries)

			SAVE_LIBS=${LIBS}
			LIBS="${QTQGL_LIBOPTS} ${OPENGL_LIBOPTS} ${X11_LIBOPTS} ${LIBS} ${VIEW_INCLUDES}"
			AC_TRY_LINK([#include <qgl.h>], [QGLWidget widget;], QT_LINKING_OK=1)
			LIBS=${SAVE_LIBS}
	
			if test "${QT_LINKING_OK+set}" != set ; then
				SAVE_LIBS=${LIBS}
				LIBS="${QT_LIBOPTS} ${OPENGL_LIBOPTS} ${X11_LIBOPTS} ${LIBS} ${VIEW_INCLUDES}"
				AC_TRY_LINK([#include <qgl.h>], [QGLWidget wid;], QT_LINKING_OK=1)
				LIBS=${SAVE_LIBS}
			else
				dnl link against qgl as well (for qt <= 2.0)
				QT_LIBOPTS="${QTQGL_LIBOPTS}"
			fi

			if test "${QT_LINKING_OK+set}" != set ; then
				SAVE_LIBS=${LIBS}
				X11_LIBOPTS="-lXrender -lfreetype ${X11_LIBOPTS}"
				LIBS="${QT_LIBOPTS} ${OPENGL_LIBOPTS} ${X11_LIBOPTS} ${LIBS} ${VIEW_INCLUDES}"
				AC_TRY_LINK([#include <qgl.h>], [QGLWidget wid;], QT_LINKING_OK=1)
				LIBS=${SAVE_LIBS}
			fi

		if test "${QT_LINKING_OK+set}" != set ; then
			AC_MSG_RESULT(no)
			AC_MSG_RESULT()
			AC_MSG_RESULT([Cannot link against libqt!])
			AC_MSG_RESULT([If QT is installed, please specify the path to the library])
			AC_MSG_RESULT([using the option --with-qt-libs=DIR or the environment variable QTDIR.])
			CF_ERROR
		else
			AC_MSG_RESULT(yes)
			
			dnl  
			dnl  identify the version of the library
			dnl
			AC_MSG_CHECKING(QT library version)
			SAVE_LIBS=${LIBS}
			LIBS="${QT_LIBOPTS} ${OPENGL_LIBOPTS} ${X11_LIBOPTS} ${LIBS}"
			if test "${OS}" = "Darwin" ; then
				DYLD_LIBRARY_PATH="${QT_LIBPATH}:${X11_LIBPATH}:${OPENGL_LIBPATH}:${DYLD_LIBRARY_PATH}"
				export DYLD_LIBRARY_PATH
				echo "DYLD_LIBRARY_PATH = ${DYLD_LIBRARY_PATH}" 1>&5
			else
				LD_LIBRARY_PATH="${QT_LIBPATH}:${X11_LIBPATH}:${OPENGL_LIBPATH}:${LD_LIBRARY_PATH}"
				export LD_LIBRARY_PATH
				echo "LD_LIBRARY_PATH = ${LD_LIBRARY_PATH}" 1>&5
			fi
			AC_TRY_RUN(
				[
					#include <stdio.h> 
					const char* qVersion();
					int main()
					{
						FILE* f = fopen("qt.version", "w");
						fprintf(f, "%s\n", qVersion());
						fclose(f);
						return 0;
					}
				], 
				QT_VERSION_OK=1,
				DUMMY=0,
				DUMMY=0
			)
			LIBS=${SAVE_LIBS}
			
			dnl
			dnl	if the program compiled and ran successfully,
			dnl extract the QT version number
			dnl
			if test "${QT_VERSION_OK+set}" != set; then
				AC_MSG_RESULT(no)
				AC_MSG_RESULT()
				AC_MSG_RESULT(The execution of a program linked against the QT)
				AC_MSG_RESULT(library failed. Please have a look at config.log)
				AC_MSG_RESULT((the last few lines) to find out what happened.)
				AC_MSG_RESULT(Perhaps you specified the wrong library or the)
				AC_MSG_RESULT(X11 libraries are in conflict with any other library.)
				AC_MSG_RESULT(You might also want to check your LD_LIBRARY_PATH.)
				CF_ERROR
			else
				QT_VERSION_STRING=`cat qt.version`
				AC_MSG_RESULT(${QT_VERSION_STRING})

				dnl
				dnl  test whether this version is the right one
				dnl  (2.x.y and at least 2.0.2
				dnl
				${RM} qt.version 2>/dev/null
				QT_MAJOR=`echo ${QT_VERSION_STRING} | ${CUT} -d. -f1`
				if test "${QT_MAJOR}" -lt 3 ; then
					AC_MSG_RESULT()
					AC_MSG_RESULT(QT version 3.x is required.)
					AC_MSG_RESULT(Please install version QT Version 3 (at least 3.0.6))
					AC_MSG_RESULT(which can be obtained from)
					AC_MSG_RESULT()
					AC_MSG_RESULT(  www.troll.no/qt)
					CF_ERROR
				fi
			fi
		fi
	fi

	dnl
	dnl	try to find the MOC (QT meta object compiler)
	dnl It is usually installed in ${QTDIR}/bin/moc
	dnl
	if test "${USE_VIEW}" = true ; then
		if test "${MOC}" = moc ; then
			if test "${QTDIR}" != "" ; then
				MOC=${QTDIR}/bin/moc
			fi
		fi

		dnl
		dnl  try to find that damned moc
		dnl
		AC_PATH_PROG(MOC,moc,moc)
		if test "${MOC}" = moc ; then
			AC_MSG_RESULT()
			AC_MSG_RESULT([Could not find the QT Meta Object Compiler (moc)!])
			AC_MSG_RESULT([You might run into trouble if you want to compile MolVIEW.])
			AC_MSG_RESULT([Please include the correct path to moc into your])
			AC_MSG_RESULT([PATH environment variable or specify the path to moc])
			AC_MSG_RESULT([using the option --with-moc=PATH to rerun configure.])
			CF_ERROR
		fi
    dnl
    dnl  Make sure the MOC we found is actually executable
    dnl
		AC_MSG_CHECKING(whether we can run moc)
    if test ! -x "${MOC}" ; then
			AC_MSG_RESULT()
			AC_MSG_RESULT([The QT Meta Object Compiler (moc) found in ])
      AC_MSG_RESULT("  ${MOC}")
      AC_MSG_RESULT([seems not to be an executable!])
			AC_MSG_RESULT([Please include the correct path to moc into your])
			AC_MSG_RESULT([PATH environment variable or specify the path to moc])
			AC_MSG_RESULT([using the option --with-moc=PATH to rerun configure.])
			CF_ERROR
		else
			AC_MSG_RESULT(yes)
			AC_MSG_CHECKING(moc version)
			MOC_VERSION=`${MOC} -v 2>&1 | ${TR} -d "()" | ${SED} "s/.* Qt //"`
			AC_MSG_RESULT(${MOC_VERSION})
			
			if test "${MOC_VERSION}" != "${QT_VERSION_STR}" ; then
				AC_MSG_RESULT()
				AC_MSG_RESULT([QT version (${QT_VERSION_STR}) is incompatible with moc version (${MOC_VERISON})!])
				AC_MSG_RESULT([Please check your QTDRI environment variable, include the correct])
				AC_MSG_RESULT([path to moc in your PATH environment variable, or specify the correct])
				AC_MSG_RESULT([path to moc using the option --with-moc=PATH to rerun configure.])
				CF_ERROR
			fi
		fi
	fi

	dnl
	dnl	try to find the UIC (QT user interface compiler)
	dnl It is usually installed in ${QTDIR}/bin/uic
	dnl
	if test "${USE_VIEW}" = true ; then
		if test "${UIC}" = uic ; then
			if test "${QTDIR}" != "" ; then
				UIC=${QTDIR}/bin/uic
			fi
		fi

		dnl
		dnl  try to find that damned uic
		dnl
		AC_PATH_PROG(UIC,uic,uic)
		if test "${UIC}" = uic ; then
			AC_MSG_RESULT()
			AC_MSG_RESULT([Could not find the QT User Interface Compiler (uic)!])
			AC_MSG_RESULT([You might run into trouble if you want to compile MolVIEW.])
			AC_MSG_RESULT([Please include the correct path to uic into your])
			AC_MSG_RESULT([PATH environment variable or specify the path to uic])
			AC_MSG_RESULT([using the option --with-uic=PATH to rerun configure.])
			CF_ERROR
		fi
    
    dnl
    dnl  Make sure the UIC we found is actually executable
    dnl
		AC_MSG_CHECKING(whether uic is executable)
    if test ! -x "${UIC}" ; then
			AC_MSG_RESULT()
			AC_MSG_RESULT([The QT User Interface Compiler (uic) found in ])
      AC_MSG_RESULT("   ${UIC}")
      AC_MSG_RESULT([seems not to be an executable!])
			AC_MSG_RESULT([Please include the correct path to uic into your])
			AC_MSG_RESULT([PATH environment variable or specify the path to uic])
			AC_MSG_RESULT([using the option --with-uic=PATH to rerun configure.])
			CF_ERROR
		else
			AC_MSG_RESULT(yes)
			AC_MSG_CHECKING(uic version)
			UIC_VERSION=`${UIC} -version 2>&1 | ${GREP} version | ${TR} -d "()" | ${SED} "s/.*version //"`
			AC_MSG_RESULT(${UIC_VERSION})
			
			if test "${UIC_VERSION}" != "${QT_VERSION_STR}" ; then
				AC_MSG_RESULT()
				AC_MSG_RESULT([QT version (${QT_VERSION_STR}) is incompatible with uic version (${UIC_VERISON})!])
				AC_MSG_RESULT([Please check your QTDIR environment variable, include the correct])
				AC_MSG_RESULT([path to uic in your PATH environment variable, or specify the correct])
				AC_MSG_RESULT([path to uic using the option --with-uic=PATH to rerun configure.])
				CF_ERROR
			fi
		fi
	fi
    

	if test "${USE_VIEW}" = "true" ; then
		AC_DEFINE(PROJECT[]_HAS_VIEW,)
		LIBVIEW="libVIEW.a"
		VIEW="VIEW"
	else
		VIEW=
	fi
])



dnl
dnl		Python extension support
dnl
AC_DEFUN(CF_PYTHON, [
	dnl
	dnl  variable subsitutions for config.mak
	dnl
	AC_SUBST(SIP)
	AC_SUBST(SIP_LIB)
	AC_SUBST(SIP_INCLUDES)
	AC_SUBST(PYTHON_INCLUDES)
	AC_SUBST(PYTHON_SUPPORT)
	AC_SUBST(PYTHON_LIBS)

	if test "${PYTHON_SUPPORT}" = true ; then
		dnl
		dnl Python support won't work without VIEW!
		dnl (at least for the moment...)
		dnl
		if test "${USE_VIEW}" = false ; then
			AC_MSG_RESULT()
			AC_MSG_RESULT(PROJECT[] Python support requires the visualization component)
			AC_MSG_RESULT(VIEW. Please reconfigure without --without-VIEW.)
			CF_ERROR
		fi
		AC_DEFINE(PROJECT[]_PYTHON_SUPPORT)

		dnl 
		dnl Find the python executable (specified via --with-python)
		dnl 
	 
		dnl
		dnl  If a complete path is specified, assume it is correct.
		dnl  Otherwise, search in the current PATH for a suitable executable.
		dnl
		if test "`basename ${PYTHON_EXE}`" = "${PYTHON_EXE}" ; then 
			AC_PATH_PROG(PYTHON_EXECUTABLE, ${PYTHON_EXE})
		else
			PYTHON_EXECUTABLE="${PYTHON_EXE}"
		fi
	 
		if test "${PYTHON_EXECUTABLE}" = "" ; then
			AC_MSG_RESULT()
			AC_MSG_RESULT([Could not find Python interpreter ]${PYTHON_EXCUTABLE})
			AC_MSG_RESULT([Please use --with-python=EXE to specify its location.])
			CF_ERROR
		fi

		AC_MSG_CHECKING(for Python interpreter)
		if test -x "${PYTHON_EXECUTABLE}" ; then
			AC_MSG_RESULT(${PYTHON_EXECUTABLE})
		else
			AC_MSG_RESULT()
			AC_MSG_RESULT([Could not find Python interpreter ]${PYTHON_EXCUTABLE})
			AC_MSG_RESULT([Please use --with-python=EXE to specify its location.])
			CF_ERROR
		fi
			
		dnl
		dnl	 Run python to retrieve some useful configuration information
		dnl	
		AC_MSG_CHECKING(for Python version)
		PYTHON_VERSION=`${PYTHON_EXECUTABLE} -c 'import sys;print sys.version' | ${SED} -n 1p | ${CUT} -d\  -f1`
		AC_MSG_RESULT(${PYTHON_VERSION})
		PYTHON_VERSION_NUMBER_1=`echo ${PYTHON_VERSION} | ${CUT} -d. -f1`
		PYTHON_VERSION_NUMBER_2=`echo ${PYTHON_VERSION} | ${CUT} -d. -f2`
		PYTHON_VERSION_NUMBER_3=`echo ${PYTHON_VERSION} | ${CUT} -d. -f3`

		dnl
		dnl	shorten the release number to Major.minor (only those are used to construct
		dnl include and lib paths)
		dnl
		PYTHON_VERSION="${PYTHON_VERSION_NUMBER_1}.${PYTHON_VERSION_NUMBER_2}"
		
		dnl
		dnl  We need at least Python 2.3
		dnl
		if test "${PYTHON_VERSION_NUMBER_1}" -le 1 -o "${PYTHON_VERSION_NUMBER_2}" -lt 2 ; then
			AC_MSG_RESULT()
			AC_MSG_RESULT([Python verison 2.3 or above required!])
			AC_MSG_RESULT([Please donwload and install Python from])
			AC_MSG_RESULT([  http://www.python.org])
			CF_ERROR
		fi
		
		AC_MSG_CHECKING(for Python installation paths)
		PYTHON_PREFIX=`${PYTHON_EXECUTABLE} -c 'import sys;print sys.prefix'`
		AC_MSG_RESULT(${PYTHON_PREFIX})


		dnl
		dnl  Python include path
		dnl
		AC_MSG_CHECKING(for Python.h)
		if test "${PYTHON_INC_PATH}" = "" ; then
			PYTHON_INC_PATH="${PYTHON_PREFIX}/include/python${PYTHON_VERSION}"
		fi
		CF_FIND_HEADER(PYTHON_INC_PATH, Python.h, ${PYTHON_INCLUDE_PATH})
		if test "${PYTHON_INC_PATH}" = "" ; then
			AC_MSG_RESULT(not found!)
			AC_MSG_RESULT()
			AC_MSG_RESULT([Please specify the path to the directory that contains])
			AC_MSG_RESULT([Python.h using the option --with-python-incl=DIR])
			AC_MSG_RESULT([or ensure that Python is installed in the correct directory])
			AC_MSG_RESULT([(sys.prefix is ${PYTHON_PREFIX})])
			CF_ERROR
		else
			AC_MSG_RESULT(${PYTHON_INC_PATH})
			PYTHON_INCLUDES="-I${PYTHON_INC_PATH}"
		fi

		dnl
		dnl	Python library path
		dnl
    dnl use framework Python instead for Darwin
    if test "${OS}" = "Darwin" ; then
 			PYTHON_LIBS="-framework Python"
		else     
			AC_MSG_CHECKING(for libpython)
			if test "${PYTHON_LIBPATH}" = "" ; then
				PYTHON_LIBPATH="${PYTHON_PREFIX}/lib/python${PYTHON_VERSION}/config/"
			fi
			PYTHON_LIBS=`${FIND} ${PYTHON_LIBPATH} -name libpython*.a 2>/dev/null`
			if test "${PYTHON_LIBS}" = "" ; then
				AC_MSG_RESULT()
				AC_MSG_RESULT(No libpython*a found in ${PYTHON_LIBPATH}. Please specify)
				AC_MSG_RESULT(the path where your Python library resides using --with-python-libs=DIR)
				AC_MSG_RESULT(or ensure that libpython is installed in the correct directory)
				AC_MSG_RESULT([(sys.prefix is ]${PYTHON_PREFIX}[)])
				CF_ERROR
			fi
			AC_MSG_RESULT(${PYTHON_LIBS})

			if test "${PYTHON_LDOPTS}" = "" ; then
				PYTHON_MAKEFILE=`${FIND} ${PYTHON_LIBPATH} -name Makefile 2>/dev/null`
				if test "${PYTHON_MAKEFILE}" = "" ; then
					AC_MSG_RESULT()
					AC_MSG_RESULT(Makefile in the Python lib/config directory not found!)
					AC_MSG_RESULT(Please specify the correct options needed to link)
					AC_MSG_RESULT(against the Python library using)
					AC_MSG_RESULT( --with-python-ldopts=OPTIONS)
					AC_MSG_RESULT([(e.g. --with-python-ldopts="-ltermcap -lm")])
					CF_ERROR
				fi
				PYTHON_LIBS="${PYTHON_LIBS} `${GREP} \^LIBS= ${PYTHON_MAKEFILE} | ${CUT} -d=  -f2-`"
				PYTHON_LIBS="${PYTHON_LIBS} `${GREP} \^BASEMODLIBS= ${PYTHON_MAKEFILE} | ${CUT} -d=  -f2-`"
				PYTHON_LIBS="${PYTHON_LIBS} `${GREP} \^LOCALMODLIBS= ${PYTHON_MAKEFILE} | ${CUT} -d=  -f2-` -lm"
				PYTHON_LIBS=`echo ${PYTHON_LIBS} | ${TR} -s " "`
			fi
			AC_MSG_RESULT(Linker options for Python library: ${PYTHON_LIBS})
		fi

		
		dnl
		dnl	 SIP executable
		dnl
		CF_MSG_PATH_PROG(SIP,sip,no)
		if test "${SIP}" = no ; then
			AC_MSG_RESULT()
			AC_MSG_RESULT(Could not find SIP executable!)
			AC_MSG_RESULT(Please specify the location of SIP using the)
			AC_MSG_RESULT( --with-sip=PATH)
			AC_MSG_RESULT(option or make sure it is in your current PATH.)
			CF_ERROR
		fi

		dnl
		dnl	 SIP executable
		dnl
		AC_MSG_CHECKING(whether ${SIP} is executable)
		if test ! -x "${SIP}" ; then
			AC_MSG_RESULT()
			AC_MSG_RESULT(Could not execute ${SIP}!)
			AC_MSG_RESULT(Please specify the location of SIP using the)
			AC_MSG_RESULT( --with-sip=PATH)
			AC_MSG_RESULT(option or make sure it is in your current PATH.)
			CF_ERROR
		fi
		AC_MSG_RESULT(yes)

		dnl
		dnl		SIP version
		dnl
		AC_MSG_CHECKING(sip version)
		SIP_VERSION=`$SIP -V`
		AC_MSG_RESULT(${SIP_VERSION})
		if test "${SIP}" = "" ; then
			AC_MSG_RESULT()
			AC_MSG_RESULT(Could not determine version number of ${SIP}!)
			AC_MSG_RESULT(Please specify the location of SIP using the)
			AC_MSG_RESULT( --with-sip=PATH)
			AC_MSG_RESULT(option or make sure it is in your current PATH.)
			CF_ERROR
		fi
		SIP_VERS_NUM=`echo ${SIP_VERSION}| ${CUT} -d\  -f1`
		SIP_VERS_MAJOR=`echo ${SIP_VERS_NUM} | ${CUT} -d. -f1`
		SIP_VERS_MINOR=`echo ${SIP_VERS_NUM} | ${CUT} -d. -f2`
		SIP_VERS_MINOR_MINOR=`echo ${SIP_VERS_NUM} | ${CUT} -d. -f3`
		if test "${SIP_VERS_MINOR_MINOR}" = "" ; then
			SIP_VERS_MINOR_MINOR="0"
		fi
		if test "${SIP_VERS_MAJOR}" -lt 4 \
				-o "${SIP_VERS_MAJOR}" = 4 -a "${SIP_VERS_MINOR}" -lt 3 \
				-o "${SIP_VERS_MAJOR}" = 4 -a "${SIP_VERS_MINOR}" = 3 -a "${SIP_VERS_MINOR_MINOR}" -lt 1; then
			AC_MSG_RESULT()
			AC_MSG_RESULT(SIP release 4.3.1 or above required.)
			AC_MSG_RESULT(Your version: ${SIP_VERSION}")
			AC_MSG_RESULT(Please upgrade or specify the location of the correct SIP using the)
			AC_MSG_RESULT( --with-sip=PATH)
			AC_MSG_RESULT()
			CF_ERROR
		fi
	
		dnl
		dnl	SIP header file (sip.h)
		dnl
		AC_MSG_CHECKING(sip header file)
		if test "${SIP_INCLUDE_PATH}" = "" ; then
			SIP_INCLUDE_PATH="${PYTHON_INC_PATH}"
		fi
		CF_FIND_HEADER(SIP_INC_PATH, sip.h, ${SIP_INCLUDE_PATH})
		if test "${SIP_INC_PATH}" = "" ; then
			AC_MSG_RESULT(not found!)
			AC_MSG_RESULT()
			AC_MSG_RESULT(Please specify the path to the directory that contains)
			AC_MSG_RESULT(sip.h using the option --with-sip-incl=DIR.)
			AC_MSG_RESULT([If you do not have that file, you should obtain SIP])
			AC_MSG_RESULT(from)
			AC_MSG_RESULT(  www.thekompany.com/projects/pykde)
			CF_ERROR
		else
			AC_MSG_RESULT(${SIP_INC_PATH})
			SIP_INCLUDES="-I${SIP_INC_PATH}"
		fi
	fi

	if test "${PYTHON_SUPPORT}" = true ; then
		PYTHON=PYTHON
	else
		PYTHON=
	fi
])


# Do all the work for Automake.                            -*- Autoconf -*-

# This macro actually does too much some checks are only needed if
# your package does certain things.  But this isn't really a big deal.

# Copyright 1996, 1997, 1998, 1999, 2000, 2001, 2002
# Free Software Foundation, Inc.

# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2, or (at your option)
# any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
# 02111-1307, USA.

# serial 8

# There are a few dirty hacks below to avoid letting `AC_PROG_CC' be
# written in clear, in which case automake, when reading aclocal.m4,
# will think it sees a *use*, and therefore will trigger all it's
# C support machinery.  Also note that it means that autoscan, seeing
# CC etc. in the Makefile, will ask for an AC_PROG_CC use...


AC_PREREQ([2.52])

# Autoconf 2.50 wants to disallow AM_ names.  We explicitly allow
# the ones we care about.
m4_pattern_allow([^AM_[A-Z]+FLAGS$])dnl

# AM_INIT_AUTOMAKE(PACKAGE, VERSION, [NO-DEFINE])
# AM_INIT_AUTOMAKE([OPTIONS])
# -----------------------------------------------
# The call with PACKAGE and VERSION arguments is the old style
# call (pre autoconf-2.50), which is being phased out.  PACKAGE
# and VERSION should now be passed to AC_INIT and removed from
# the call to AM_INIT_AUTOMAKE.
# We support both call styles for the transition.  After
# the next Automake release, Autoconf can make the AC_INIT
# arguments mandatory, and then we can depend on a new Autoconf
# release and drop the old call support.
AC_DEFUN([AM_INIT_AUTOMAKE],
[AC_REQUIRE([AM_SET_CURRENT_AUTOMAKE_VERSION])dnl
 AC_REQUIRE([AC_PROG_INSTALL])dnl
# test to see if srcdir already configured
if test "`cd $srcdir && pwd`" != "`pwd`" &&
   test -f $srcdir/config.status; then
  AC_MSG_ERROR([source directory already configured; run "make distclean" there first])
fi

# Define the identity of the package.
dnl Distinguish between old-style and new-style calls.
m4_ifval([$2],
[m4_ifval([$3], [_AM_SET_OPTION([no-define])])dnl
 AC_SUBST([PACKAGE], [$1])dnl
 AC_SUBST([VERSION], [$2])],
[_AM_SET_OPTIONS([$1])dnl
 AC_SUBST([PACKAGE], [AC_PACKAGE_TARNAME])dnl
 AC_SUBST([VERSION], [AC_PACKAGE_VERSION])])dnl

_AM_IF_OPTION([no-define],,
[AC_DEFINE_UNQUOTED(PACKAGE, "$PACKAGE", [Name of package])
 AC_DEFINE_UNQUOTED(VERSION, "$VERSION", [Version number of package])])dnl

# Some tools Automake needs.
AC_REQUIRE([AM_SANITY_CHECK])dnl
AC_REQUIRE([AC_ARG_PROGRAM])dnl
AM_MISSING_PROG(ACLOCAL, aclocal-${am__api_version})
AM_MISSING_PROG(AUTOCONF, autoconf)
AM_MISSING_PROG(AUTOMAKE, automake-${am__api_version})
AM_MISSING_PROG(AUTOHEADER, autoheader)
AM_MISSING_PROG(MAKEINFO, makeinfo)
AM_MISSING_PROG(AMTAR, tar)
AM_PROG_INSTALL_SH
AM_PROG_INSTALL_STRIP
# We need awk for the "check" target.  The system "awk" is bad on
# some platforms.
AC_REQUIRE([AC_PROG_AWK])dnl
AC_REQUIRE([AC_PROG_MAKE_SET])dnl

_AM_IF_OPTION([no-dependencies],,
[AC_PROVIDE_IFELSE([AC_PROG_][CC],
                  [_AM_DEPENDENCIES(CC)],
                  [define([AC_PROG_][CC],
                          defn([AC_PROG_][CC])[_AM_DEPENDENCIES(CC)])])dnl
AC_PROVIDE_IFELSE([AC_PROG_][CXX],
                  [_AM_DEPENDENCIES(CXX)],
                  [define([AC_PROG_][CXX],
                          defn([AC_PROG_][CXX])[_AM_DEPENDENCIES(CXX)])])dnl
])
])

# Copyright 2002  Free Software Foundation, Inc.

# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2, or (at your option)
# any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA

# AM_AUTOMAKE_VERSION(VERSION)
# ----------------------------
# Automake X.Y traces this macro to ensure aclocal.m4 has been
# generated from the m4 files accompanying Automake X.Y.
AC_DEFUN([AM_AUTOMAKE_VERSION],[am__api_version="1.6"])

# AM_SET_CURRENT_AUTOMAKE_VERSION
# -------------------------------
# Call AM_AUTOMAKE_VERSION so it can be traced.
# This function is AC_REQUIREd by AC_INIT_AUTOMAKE.
AC_DEFUN([AM_SET_CURRENT_AUTOMAKE_VERSION],
	 [AM_AUTOMAKE_VERSION([1.6.3])])

# Helper functions for option handling.                    -*- Autoconf -*-

# Copyright 2001, 2002  Free Software Foundation, Inc.

# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2, or (at your option)
# any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
# 02111-1307, USA.

# serial 2

# _AM_MANGLE_OPTION(NAME)
# -----------------------
AC_DEFUN([_AM_MANGLE_OPTION],
[[_AM_OPTION_]m4_bpatsubst($1, [[^a-zA-Z0-9_]], [_])])

# _AM_SET_OPTION(NAME)
# ------------------------------
# Set option NAME.  Presently that only means defining a flag for this option.
AC_DEFUN([_AM_SET_OPTION],
[m4_define(_AM_MANGLE_OPTION([$1]), 1)])

# _AM_SET_OPTIONS(OPTIONS)
# ----------------------------------
# OPTIONS is a space-separated list of Automake options.
AC_DEFUN([_AM_SET_OPTIONS],
[AC_FOREACH([_AM_Option], [$1], [_AM_SET_OPTION(_AM_Option)])])

# _AM_IF_OPTION(OPTION, IF-SET, [IF-NOT-SET])
# -------------------------------------------
# Execute IF-SET if OPTION is set, IF-NOT-SET otherwise.
AC_DEFUN([_AM_IF_OPTION],
[m4_ifset(_AM_MANGLE_OPTION([$1]), [$2], [$3])])

#
# Check to make sure that the build environment is sane.
#

# Copyright 1996, 1997, 2000, 2001 Free Software Foundation, Inc.

# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2, or (at your option)
# any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
# 02111-1307, USA.

# serial 3

# AM_SANITY_CHECK
# ---------------
AC_DEFUN([AM_SANITY_CHECK],
[AC_MSG_CHECKING([whether build environment is sane])
# Just in case
sleep 1
echo timestamp > conftest.file
# Do `set' in a subshell so we don't clobber the current shell's
# arguments.  Must try -L first in case configure is actually a
# symlink; some systems play weird games with the mod time of symlinks
# (eg FreeBSD returns the mod time of the symlink's containing
# directory).
if (
   set X `ls -Lt $srcdir/configure conftest.file 2> /dev/null`
   if test "$[*]" = "X"; then
      # -L didn't work.
      set X `ls -t $srcdir/configure conftest.file`
   fi
   rm -f conftest.file
   if test "$[*]" != "X $srcdir/configure conftest.file" \
      && test "$[*]" != "X conftest.file $srcdir/configure"; then

      # If neither matched, then we have a broken ls.  This can happen
      # if, for instance, CONFIG_SHELL is bash and it inherits a
      # broken ls alias from the environment.  This has actually
      # happened.  Such a system could not be considered "sane".
      AC_MSG_ERROR([ls -t appears to fail.  Make sure there is not a broken
alias in your environment])
   fi

   test "$[2]" = conftest.file
   )
then
   # Ok.
   :
else
   AC_MSG_ERROR([newly created file is older than distributed files!
Check your system clock])
fi
AC_MSG_RESULT(yes)])

#  -*- Autoconf -*-


# Copyright 1997, 1999, 2000, 2001 Free Software Foundation, Inc.

# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2, or (at your option)
# any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
# 02111-1307, USA.

# serial 3

# AM_MISSING_PROG(NAME, PROGRAM)
# ------------------------------
AC_DEFUN([AM_MISSING_PROG],
[AC_REQUIRE([AM_MISSING_HAS_RUN])
$1=${$1-"${am_missing_run}$2"}
AC_SUBST($1)])


# AM_MISSING_HAS_RUN
# ------------------
# Define MISSING if not defined so far and test if it supports --run.
# If it does, set am_missing_run to use it, otherwise, to nothing.
AC_DEFUN([AM_MISSING_HAS_RUN],
[AC_REQUIRE([AM_AUX_DIR_EXPAND])dnl
test x"${MISSING+set}" = xset || MISSING="\${SHELL} $am_aux_dir/missing"
# Use eval to expand $SHELL
if eval "$MISSING --run true"; then
  am_missing_run="$MISSING --run "
else
  am_missing_run=
  AC_MSG_WARN([`missing' script is too old or missing])
fi
])

# AM_AUX_DIR_EXPAND

# Copyright 2001 Free Software Foundation, Inc.

# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2, or (at your option)
# any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
# 02111-1307, USA.

# For projects using AC_CONFIG_AUX_DIR([foo]), Autoconf sets
# $ac_aux_dir to `$srcdir/foo'.  In other projects, it is set to
# `$srcdir', `$srcdir/..', or `$srcdir/../..'.
#
# Of course, Automake must honor this variable whenever it calls a
# tool from the auxiliary directory.  The problem is that $srcdir (and
# therefore $ac_aux_dir as well) can be either absolute or relative,
# depending on how configure is run.  This is pretty annoying, since
# it makes $ac_aux_dir quite unusable in subdirectories: in the top
# source directory, any form will work fine, but in subdirectories a
# relative path needs to be adjusted first.
#
# $ac_aux_dir/missing
#    fails when called from a subdirectory if $ac_aux_dir is relative
# $top_srcdir/$ac_aux_dir/missing
#    fails if $ac_aux_dir is absolute,
#    fails when called from a subdirectory in a VPATH build with
#          a relative $ac_aux_dir
#
# The reason of the latter failure is that $top_srcdir and $ac_aux_dir
# are both prefixed by $srcdir.  In an in-source build this is usually
# harmless because $srcdir is `.', but things will broke when you
# start a VPATH build or use an absolute $srcdir.
#
# So we could use something similar to $top_srcdir/$ac_aux_dir/missing,
# iff we strip the leading $srcdir from $ac_aux_dir.  That would be:
#   am_aux_dir='\$(top_srcdir)/'`expr "$ac_aux_dir" : "$srcdir//*\(.*\)"`
# and then we would define $MISSING as
#   MISSING="\${SHELL} $am_aux_dir/missing"
# This will work as long as MISSING is not called from configure, because
# unfortunately $(top_srcdir) has no meaning in configure.
# However there are other variables, like CC, which are often used in
# configure, and could therefore not use this "fixed" $ac_aux_dir.
#
# Another solution, used here, is to always expand $ac_aux_dir to an
# absolute PATH.  The drawback is that using absolute paths prevent a
# configured tree to be moved without reconfiguration.

# Rely on autoconf to set up CDPATH properly.
AC_PREREQ([2.50])

AC_DEFUN([AM_AUX_DIR_EXPAND], [
# expand $ac_aux_dir to an absolute path
am_aux_dir=`cd $ac_aux_dir && pwd`
])

# AM_PROG_INSTALL_SH
# ------------------
# Define $install_sh.

# Copyright 2001 Free Software Foundation, Inc.

# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2, or (at your option)
# any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
# 02111-1307, USA.

AC_DEFUN([AM_PROG_INSTALL_SH],
[AC_REQUIRE([AM_AUX_DIR_EXPAND])dnl
install_sh=${install_sh-"$am_aux_dir/install-sh"}
AC_SUBST(install_sh)])

# AM_PROG_INSTALL_STRIP

# Copyright 2001 Free Software Foundation, Inc.

# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2, or (at your option)
# any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
# 02111-1307, USA.

# One issue with vendor `install' (even GNU) is that you can't
# specify the program used to strip binaries.  This is especially
# annoying in cross-compiling environments, where the build's strip
# is unlikely to handle the host's binaries.
# Fortunately install-sh will honor a STRIPPROG variable, so we
# always use install-sh in `make install-strip', and initialize
# STRIPPROG with the value of the STRIP variable (set by the user).
AC_DEFUN([AM_PROG_INSTALL_STRIP],
[AC_REQUIRE([AM_PROG_INSTALL_SH])dnl
# Installed binaries are usually stripped using `strip' when the user
# run `make install-strip'.  However `strip' might not be the right
# tool to use in cross-compilation environments, therefore Automake
# will honor the `STRIP' environment variable to overrule this program.
dnl Don't test for $cross_compiling = yes, because it might be `maybe'.
if test "$cross_compiling" != no; then
  AC_CHECK_TOOL([STRIP], [strip], :)
fi
INSTALL_STRIP_PROGRAM="\${SHELL} \$(install_sh) -c -s"
AC_SUBST([INSTALL_STRIP_PROGRAM])])

# serial 4						-*- Autoconf -*-

# Copyright 1999, 2000, 2001 Free Software Foundation, Inc.

# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2, or (at your option)
# any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
# 02111-1307, USA.


# There are a few dirty hacks below to avoid letting `AC_PROG_CC' be
# written in clear, in which case automake, when reading aclocal.m4,
# will think it sees a *use*, and therefore will trigger all it's
# C support machinery.  Also note that it means that autoscan, seeing
# CC etc. in the Makefile, will ask for an AC_PROG_CC use...



# _AM_DEPENDENCIES(NAME)
# ----------------------
# See how the compiler implements dependency checking.
# NAME is "CC", "CXX", "GCJ", or "OBJC".
# We try a few techniques and use that to set a single cache variable.
#
# We don't AC_REQUIRE the corresponding AC_PROG_CC since the latter was
# modified to invoke _AM_DEPENDENCIES(CC); we would have a circular
# dependency, and given that the user is not expected to run this macro,
# just rely on AC_PROG_CC.
AC_DEFUN([_AM_DEPENDENCIES],
[AC_REQUIRE([AM_SET_DEPDIR])dnl
AC_REQUIRE([AM_OUTPUT_DEPENDENCY_COMMANDS])dnl
AC_REQUIRE([AM_MAKE_INCLUDE])dnl
AC_REQUIRE([AM_DEP_TRACK])dnl

ifelse([$1], CC,   [depcc="$CC"   am_compiler_list=],
       [$1], CXX,  [depcc="$CXX"  am_compiler_list=],
       [$1], OBJC, [depcc="$OBJC" am_compiler_list='gcc3 gcc'],
       [$1], GCJ,  [depcc="$GCJ"  am_compiler_list='gcc3 gcc'],
                   [depcc="$$1"   am_compiler_list=])

AC_CACHE_CHECK([dependency style of $depcc],
               [am_cv_$1_dependencies_compiler_type],
[if test -z "$AMDEP_TRUE" && test -f "$am_depcomp"; then
  # We make a subdir and do the tests there.  Otherwise we can end up
  # making bogus files that we don't know about and never remove.  For
  # instance it was reported that on HP-UX the gcc test will end up
  # making a dummy file named `D' -- because `-MD' means `put the output
  # in D'.
  mkdir conftest.dir
  # Copy depcomp to subdir because otherwise we won't find it if we're
  # using a relative directory.
  cp "$am_depcomp" conftest.dir
  cd conftest.dir

  am_cv_$1_dependencies_compiler_type=none
  if test "$am_compiler_list" = ""; then
     am_compiler_list=`sed -n ['s/^#*\([a-zA-Z0-9]*\))$/\1/p'] < ./depcomp`
  fi
  for depmode in $am_compiler_list; do
    # We need to recreate these files for each test, as the compiler may
    # overwrite some of them when testing with obscure command lines.
    # This happens at least with the AIX C compiler.
    echo '#include "conftest.h"' > conftest.c
    echo 'int i;' > conftest.h
    echo "${am__include} ${am__quote}conftest.Po${am__quote}" > confmf

    case $depmode in
    nosideeffect)
      # after this tag, mechanisms are not by side-effect, so they'll
      # only be used when explicitly requested
      if test "x$enable_dependency_tracking" = xyes; then
	continue
      else
	break
      fi
      ;;
    none) break ;;
    esac
    # We check with `-c' and `-o' for the sake of the "dashmstdout"
    # mode.  It turns out that the SunPro C++ compiler does not properly
    # handle `-M -o', and we need to detect this.
    if depmode=$depmode \
       source=conftest.c object=conftest.o \
       depfile=conftest.Po tmpdepfile=conftest.TPo \
       $SHELL ./depcomp $depcc -c conftest.c -o conftest.o >/dev/null 2>&1 &&
       grep conftest.h conftest.Po > /dev/null 2>&1 &&
       ${MAKE-make} -s -f confmf > /dev/null 2>&1; then
      am_cv_$1_dependencies_compiler_type=$depmode
      break
    fi
  done

  cd ..
  rm -rf conftest.dir
else
  am_cv_$1_dependencies_compiler_type=none
fi
])
AC_SUBST([$1DEPMODE], [depmode=$am_cv_$1_dependencies_compiler_type])
])


# AM_SET_DEPDIR
# -------------
# Choose a directory name for dependency files.
# This macro is AC_REQUIREd in _AM_DEPENDENCIES
AC_DEFUN([AM_SET_DEPDIR],
[rm -f .deps 2>/dev/null
mkdir .deps 2>/dev/null
if test -d .deps; then
  DEPDIR=.deps
else
  # MS-DOS does not allow filenames that begin with a dot.
  DEPDIR=_deps
fi
rmdir .deps 2>/dev/null
AC_SUBST([DEPDIR])
])


# AM_DEP_TRACK
# ------------
AC_DEFUN([AM_DEP_TRACK],
[AC_ARG_ENABLE(dependency-tracking,
[  --disable-dependency-tracking Speeds up one-time builds
  --enable-dependency-tracking  Do not reject slow dependency extractors])
if test "x$enable_dependency_tracking" != xno; then
  am_depcomp="$ac_aux_dir/depcomp"
  AMDEPBACKSLASH='\'
fi
AM_CONDITIONAL([AMDEP], [test "x$enable_dependency_tracking" != xno])
AC_SUBST([AMDEPBACKSLASH])
])

# Generate code to set up dependency tracking.   -*- Autoconf -*-

# Copyright 1999, 2000, 2001, 2002 Free Software Foundation, Inc.

# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2, or (at your option)
# any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
# 02111-1307, USA.

#serial 2

# _AM_OUTPUT_DEPENDENCY_COMMANDS
# ------------------------------
AC_DEFUN([_AM_OUTPUT_DEPENDENCY_COMMANDS],
[for mf in $CONFIG_FILES; do
  # Strip MF so we end up with the name of the file.
  mf=`echo "$mf" | sed -e 's/:.*$//'`
  # Check whether this is an Automake generated Makefile or not.
  # We used to match only the files named `Makefile.in', but
  # some people rename them; so instead we look at the file content.
  # Grep'ing the first line is not enough: some people post-process
  # each Makefile.in and add a new line on top of each file to say so.
  # So let's grep whole file.
  if grep '^#.*generated by automake' $mf > /dev/null 2>&1; then
    dirpart=`AS_DIRNAME("$mf")`
  else
    continue
  fi
  grep '^DEP_FILES *= *[[^ @%:@]]' < "$mf" > /dev/null || continue
  # Extract the definition of DEP_FILES from the Makefile without
  # running `make'.
  DEPDIR=`sed -n -e '/^DEPDIR = / s///p' < "$mf"`
  test -z "$DEPDIR" && continue
  # When using ansi2knr, U may be empty or an underscore; expand it
  U=`sed -n -e '/^U = / s///p' < "$mf"`
  test -d "$dirpart/$DEPDIR" || mkdir "$dirpart/$DEPDIR"
  # We invoke sed twice because it is the simplest approach to
  # changing $(DEPDIR) to its actual value in the expansion.
  for file in `sed -n -e '
    /^DEP_FILES = .*\\\\$/ {
      s/^DEP_FILES = //
      :loop
	s/\\\\$//
	p
	n
	/\\\\$/ b loop
      p
    }
    /^DEP_FILES = / s/^DEP_FILES = //p' < "$mf" | \
       sed -e 's/\$(DEPDIR)/'"$DEPDIR"'/g' -e 's/\$U/'"$U"'/g'`; do
    # Make sure the directory exists.
    test -f "$dirpart/$file" && continue
    fdir=`AS_DIRNAME(["$file"])`
    AS_MKDIR_P([$dirpart/$fdir])
    # echo "creating $dirpart/$file"
    echo '# dummy' > "$dirpart/$file"
  done
done
])# _AM_OUTPUT_DEPENDENCY_COMMANDS


# AM_OUTPUT_DEPENDENCY_COMMANDS
# -----------------------------
# This macro should only be invoked once -- use via AC_REQUIRE.
#
# This code is only required when automatic dependency tracking
# is enabled.  FIXME.  This creates each `.P' file that we will
# need in order to bootstrap the dependency handling code.
AC_DEFUN([AM_OUTPUT_DEPENDENCY_COMMANDS],
[AC_CONFIG_COMMANDS([depfiles],
     [test x"$AMDEP_TRUE" != x"" || _AM_OUTPUT_DEPENDENCY_COMMANDS],
     [AMDEP_TRUE="$AMDEP_TRUE" ac_aux_dir="$ac_aux_dir"])
])

# Copyright 2001 Free Software Foundation, Inc.             -*- Autoconf -*-

# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2, or (at your option)
# any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
# 02111-1307, USA.

# serial 2

# AM_MAKE_INCLUDE()
# -----------------
# Check to see how make treats includes.
AC_DEFUN([AM_MAKE_INCLUDE],
[am_make=${MAKE-make}
cat > confinc << 'END'
doit:
	@echo done
END
# If we don't find an include directive, just comment out the code.
AC_MSG_CHECKING([for style of include used by $am_make])
am__include="#"
am__quote=
_am_result=none
# First try GNU make style include.
echo "include confinc" > confmf
# We grep out `Entering directory' and `Leaving directory'
# messages which can occur if `w' ends up in MAKEFLAGS.
# In particular we don't look at `^make:' because GNU make might
# be invoked under some other name (usually "gmake"), in which
# case it prints its new name instead of `make'.
if test "`$am_make -s -f confmf 2> /dev/null | fgrep -v 'ing directory'`" = "done"; then
   am__include=include
   am__quote=
   _am_result=GNU
fi
# Now try BSD make style include.
if test "$am__include" = "#"; then
   echo '.include "confinc"' > confmf
   if test "`$am_make -s -f confmf 2> /dev/null`" = "done"; then
      am__include=.include
      am__quote="\""
      _am_result=BSD
   fi
fi
AC_SUBST(am__include)
AC_SUBST(am__quote)
AC_MSG_RESULT($_am_result)
rm -f confinc confmf
])

# AM_CONDITIONAL                                              -*- Autoconf -*-

# Copyright 1997, 2000, 2001 Free Software Foundation, Inc.

# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2, or (at your option)
# any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
# 02111-1307, USA.

# serial 5

AC_PREREQ(2.52)

# AM_CONDITIONAL(NAME, SHELL-CONDITION)
# -------------------------------------
# Define a conditional.
AC_DEFUN([AM_CONDITIONAL],
[ifelse([$1], [TRUE],  [AC_FATAL([$0: invalid condition: $1])],
        [$1], [FALSE], [AC_FATAL([$0: invalid condition: $1])])dnl
AC_SUBST([$1_TRUE])
AC_SUBST([$1_FALSE])
if $2; then
  $1_TRUE=
  $1_FALSE='#'
else
  $1_TRUE='#'
  $1_FALSE=
fi
AC_CONFIG_COMMANDS_PRE(
[if test -z "${$1_TRUE}" && test -z "${$1_FALSE}"; then
  AC_MSG_ERROR([conditional \"$1\" was never defined.
Usually this means the macro was only invoked conditionally.])
fi])])

# Like AC_CONFIG_HEADER, but automatically create stamp file. -*- Autoconf -*-

# Copyright 1996, 1997, 2000, 2001 Free Software Foundation, Inc.

# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2, or (at your option)
# any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
# 02111-1307, USA.

AC_PREREQ([2.52])

# serial 6

# When config.status generates a header, we must update the stamp-h file.
# This file resides in the same directory as the config header
# that is generated.  We must strip everything past the first ":",
# and everything past the last "/".

# _AM_DIRNAME(PATH)
# -----------------
# Like AS_DIRNAME, only do it during macro expansion
AC_DEFUN([_AM_DIRNAME],
       [m4_if(regexp([$1], [^.*[^/]//*[^/][^/]*/*$]), -1,
	      m4_if(regexp([$1], [^//\([^/]\|$\)]), -1,
		    m4_if(regexp([$1], [^/.*]), -1,
			  [.],
			  patsubst([$1], [^\(/\).*], [\1])),
		    patsubst([$1], [^\(//\)\([^/].*\|$\)], [\1])),
	      patsubst([$1], [^\(.*[^/]\)//*[^/][^/]*/*$], [\1]))[]dnl
])# _AM_DIRNAME


# The stamp files are numbered to have different names.
# We could number them on a directory basis, but that's additional
# complications, let's have a unique counter.
m4_define([_AM_STAMP_Count], [0])


# _AM_STAMP(HEADER)
# -----------------
# The name of the stamp file for HEADER.
AC_DEFUN([_AM_STAMP],
[m4_define([_AM_STAMP_Count], m4_incr(_AM_STAMP_Count))dnl
AS_ESCAPE(_AM_DIRNAME(patsubst([$1],
                               [:.*])))/stamp-h[]_AM_STAMP_Count])


# _AM_CONFIG_HEADER(HEADER[:SOURCES], COMMANDS, INIT-COMMANDS)
# ------------------------------------------------------------
# We used to try to get a real timestamp in stamp-h.  But the fear is that
# that will cause unnecessary cvs conflicts.
AC_DEFUN([_AM_CONFIG_HEADER],
[# Add the stamp file to the list of files AC keeps track of,
# along with our hook.
AC_CONFIG_HEADERS([$1],
                  [# update the timestamp
echo 'timestamp for $1' >"_AM_STAMP([$1])"
$2],
                  [$3])
])# _AM_CONFIG_HEADER


# AM_CONFIG_HEADER(HEADER[:SOURCES]..., COMMANDS, INIT-COMMANDS)
# --------------------------------------------------------------
AC_DEFUN([AM_CONFIG_HEADER],
[AC_FOREACH([_AM_File], [$1], [_AM_CONFIG_HEADER(_AM_File, [$2], [$3])])
])# AM_CONFIG_HEADER


# Copyright 1998, 1999, 2000, 2001, 2002 Free Software Foundation, Inc.

# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2, or (at your option)
# any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
# 02111-1307, USA.

# serial 3

AC_PREREQ(2.50)

# AM_PROG_LEX
# -----------
# Autoconf leaves LEX=: if lex or flex can't be found.  Change that to a
# "missing" invocation, for better error output.
AC_DEFUN([AM_PROG_LEX],
[AC_REQUIRE([AM_MISSING_HAS_RUN])dnl
AC_REQUIRE([AC_PROG_LEX])dnl
if test "$LEX" = :; then
  LEX=${am_missing_run}flex
fi])


AC_DEFUN(CF_CHECK_MULTI_BUILD,[
	if test "${MULTI_BUILD}" = "true" ; then
		AC_MSG_CHECKING(multi-platform build)
		AC_MSG_RESULT(enabled)

		dnl   add the binary format to the list of supported binary formats
		dnl   held in config/binary_formats. Avoid double entries
		dnl
		if test "${MULTI_BUILD}" = "true" ; then
			touch ${BINFORMAT_FILE}
			if test "`${GREP} \^${BINFMT}\\$ ${BINFORMAT_FILE}`" = "" ; then
				echo ${BINFMT} >> ${BINFORMAT_FILE}
			fi
		fi

		dnl
		dnl   create the global config.h (the one including the platform specific
		dnl   config.h.${BINFMT})
		dnl
		${CAT} config/config.h.header | ${SED} 1,2d > config.h

		dnl
		dnl add an error line to catch all compilations without -DBMFT=
		dnl (this is usually a problem with a missing "include config.mak" in the makefile.
		dnl
		echo "#ifndef BFMT" >> config.h
		echo ["# error] PROJECT [was configured in MULTI BUILD mode! Please specify -DBMFT!" >> config.h]
		echo "#endif" >> config.h
		echo "" >> config.h

		LINES=`cat config/binary_formats | wc -l`
		i=1
		while test $i -le $LINES ; do
			BFMT=`cat ${BINFORMAT_FILE} | ${SED} -n ${i}p`
			echo "#if ( BFMT == $i )" >> config.h
			echo ["# include <]PROJECT[/CONFIG/config.h.${BFMT}>" >> config.h]
			echo "#endif" >> config.h
			echo " " >> config.h
			i=`expr $i + 1`
		done
		${CAT} config/config.h.footer | ${SED} 1,2d >> config.h
		${MKDIR} ${[]PROJECT[]_PATH}/include/[]PROJECT[]/CONFIG 2>/dev/null
		if test -f ${[]PROJECT[]_PATH}/include/[]PROJECT[]/CONFIG/config.h ; then
			if test "`${DIFF} ${[]PROJECT[]_PATH}/include/[]PROJECT[]/CONFIG/config.h config.h`" != "" ; then
				${RM} ${[]PROJECT[]_PATH}/include/PROJECT/CONFIG/config.h
				${MV} config.h  ${[]PROJECT[]_PATH}/include/[]PROJECT[]/CONFIG/config.h
			else
				${RM} config.h
			fi
		else
			${MV} config.h  ${PROJECT[]_PATH}/include/PROJECT[]/CONFIG/config.h
		fi

		dnl   define the string to substitute in common.mak
		BINFMT_PATH="/${BINFMT}"
		BINFMT_INDEX="-DBFMT="`${GREP} -n ${BINFMT} ${BINFORMAT_FILE} | ${CUT} -d: -f1 | ${TAIL} -1`
	else
		BINFMT_INDEX=""
		BINFMT_PATH=""
	fi
])

AC_DEFUN(CF_MULTI_BUILD_SHADOW, [
	if test "${MULTI_BUILD}" = "true" ; then
		AC_MSG_RESULT(creating shadow directories...)
		config/shadowsource.sh `pwd`"/${BINFMT}" `pwd` "${SUBDIRS} TEST BENCHMARKS EXAMPLES TUTORIAL APPLICATIONS"
		${RM} -fr `pwd`/${BINFMT}/TEST/data 2>/dev/null
		${RM} -fr `pwd`/${BINFMT}/BENCHMARKS/data 2>/dev/null
		${LN} -s `pwd`/TEST/data `pwd`/${BINFMT}/TEST 2>/dev/null
		${LN} -s `pwd`/TEST/runtests `pwd`/${BINFMT}/TEST 2>/dev/null
		${LN} -s `pwd`/BENCHMARKS/data `pwd`/${BINFMT}/BENCHMARKS 2>/dev/null
		${LN} -s `pwd`/BENCHMARKS/runbenchmarks `pwd`/${BINFMT}/BENCHMARKS 2>/dev/null

		${CP} config/Makefile.multiplatform Makefile
	fi
])

AC_DEFUN(CF_MOVE_CONFIG_FILES, [
	if test "${MULTI_BUILD}" = "true" ; then
		${MV} Makefile.tmp ${BINFMT}/Makefile
		${MV} common.mak.tmp ${BINFMT}/common.mak
		${MV} config.mak.tmp ${BINFMT}/config.mak
	  mkdir ${PROJECT[]_PATH}/include/PROJECT[]/CONFIG 2>/dev/null
	  ${MV} -f config.h $PROJECT[]_PATH/include/PROJECT[]/CONFIG/config.h.${BINFMT}
	else
		${MV} Makefile.tmp Makefile
		${MV} common.mak.tmp common.mak
		${MV} config.mak.tmp config.mak

		dnl
		dnl move that damned file only if it differs from the previous
		dnl version. Otherwise we have to rebuild _everything_ after each configure
		dnl
		if test -f $PROJECT[]_PATH/include/PROJECT[]/CONFIG/config.h ; then
			if test "`${DIFF} config.h $PROJECT[]_PATH/include/PROJECT[]/CONFIG/config.h`" != "" ; then
				${MV} -f config.h $PROJECT[]_PATH/include/PROJECT[]/CONFIG/config.h
			fi
		else
			dnl
			dnl  create the directory PROJECT[]/include/CONFIG
			dnl  and move config.h to that directory
			dnl
			mkdir ${PROJECT[]_PATH}/include/PROJECT[]/CONFIG 2>/dev/null
			${MV} -f config.h $PROJECT[]_PATH/include/PROJECT[]/CONFIG/config.h
		fi
	fi
])

AC_DEFUN(CF_CLEAR_DEP_FILES, [
	dnl
	dnl   make sure the dependencies and object lists are (re)built
	dnl
	if test "${MULTI_BUILD}" = "true" ; then
		${RM}  ${BINFMT}/.Dependencies 2>/dev/null
		${RM}  ${BINFMT}/lib*.objects 2>/dev/null
	else
		${RM}  .Dependencies 2>/dev/null
		${RM}  lib*.objects 2>/dev/null
	fi
])


AC_DEFUN(CF_VALGRIND, [
	dnl	
	dnl	Check for the valgrind application (a memory leak tester).
	dnl Valgrind can be used to identify leaks from the test programs
 	dnl	(target valgrind in PROJECT[]/source/TEST).
	dnl
	AC_PATH_PROG(VALGRIND, valgrind, valgrind)
	AC_SUBST(VALGRIND, $VALGRIND)
	if test "${VALGRIND}" != "valgrind" ; then
		AC_MSG_CHECKING(valgrind version)

		VALGRIND_VERSION=`${VALGRIND} --version | tr -d "a-zA-Z-_" 2>&1`
		VALGRIND_VERS_NUM=`echo ${VALGRIND_VERSION}| ${CUT} -d\  -f1`
		VALGRIND_VERS_MAJOR=`echo ${VALGRIND_VERS_NUM} | ${CUT} -d. -f1`
		VALGRIND_VERS_MINOR=`echo ${VALGRIND_VERS_NUM} | ${CUT} -d. -f2`
		VALGRIND_VERS_MINOR_MINOR=`echo ${VALGRIND_VERS_NUM} | ${CUT} -d. -f3`
		AC_MSG_RESULT(${VALGRIND_VERSION} (${VALGRIND_VERS_MAJOR}.${VALGRIND_VERS_MINOR}))	
		

		if test "${VALGRIND_VERS_MAJOR}" = "2" -a "${VALGRIND_VERS_MINOR}" -gt "0" ; then 
			AC_SUBST(VALGRIND_OPTS, "--tool=memcheck  --num-callers=20 --show-below-main=yes -v --leak-check=yes --leak-resolution=high --show-reachable=yes")
		else
			AC_SUBST(VALGRIND_OPTS, "-v --leak-check=yes --leak-resolution=high")
		fi
	else
		AC_SUBST(VALGRIND_OPTS, "-v --leak-check=yes --leak-resolution=high")
	fi
])
