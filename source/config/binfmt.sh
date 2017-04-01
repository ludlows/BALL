#!/bin/sh

############################################################
# OS and architecture detection script                     #
############################################################
# $Id: binfmt.sh,v 1.10 2000/03/27 20:05:17 oliver Exp $    #
############################################################

UNAME=uname
CUT=cut
SED=sed
GREP=egrep
TAIL=tail
EXE_DIR=`dirname $0`
BINFORMAT_FILE=${EXE_DIR}/binary_formats

OS=`${UNAME} -s`
OSREV=`${UNAME} -r`
OSMAJOR=`echo $OSREV|${CUT} -d"." -f1`

#		default...
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
	if test "${PROCESSOR}" = alpha ; then
		ARCHITECTURE=alpha
		BINFMT=Linux-alpha
	fi

	if test "${ARCHITECTURE}" = "unknown" ; then
		echo "OS: ${OS} / hardware: ${PROCESSOR}" >&2
		echo "Sorry - this architecture is currently not supported..." >&2
		exit
	fi
fi

if test ${OS} = IRIX64 ; then
	OS=IRIX
fi

if test $OS = IRIX ; then
	BINFMT=IRIX-${OSREV}
fi

if test "$OS" != Linux && test "$OS" != Solaris && test "$OS" != IRIX && test "$OS" != OSF1 && test "$OS" != FreeBSD ; then
	echo "Sorry - your OS is currently not supported..." >&2
	exit
fi

#
# 	create OS defines in config.h:
#
if test "$OS" = Linux ; then
	LINUX=LINUX
fi
if test "$OS" = Solaris ; then
	SOLARIS=SOLARIS
fi
if test "$OS" = IRIX ; then
	IRIX=IRIX
fi
if test "$OS" = OSF1 ; then
	OSF1=OSF1
fi
if test "$OS" = FreeBSD ; then
	FreeBSD=FreeBSD
fi

#
#		create ARCHITECTURE defines
#
if test "$ARCHITECTURE" = sparc ; then
	SPARC=SPARC
fi
if test "$ARCHITECTURE" = i386 ; then
	I386=I386
fi
if test "$ARCHITECTURE" = mips ; then
	MIPS=MIPS
fi
if test "$ARCHITECTURE" = alpha ; then
	ALPHA=ALPHA
fi

if test ! -f "${BINFORMAT_FILE}" ; then
	echo "cannot open file ${BINFORMAT_FILE}" >&2	
fi

FORMATS=`${GREP} "^${BINFMT}-" ${BINFORMAT_FILE}`
if test "${FORMATS}" = "" ; then
	FORMAT="undefined"
else
	if test "${FORMATS}" != "`echo ${FORMATS} | cut -d\  -f1`" ; then
		if test "${COMPILER_NAME}" != "" ; then
			COMPILER_NAME="`echo ${COMPILER_NAME} | ${SED} s/\\\\+/\\\\\\\\+/g`"
			FORMATS=`${GREP} "^${BINFMT}-" ${BINFORMAT_FILE} | ${GREP} "[-]${COMPILER_NAME}\$"`
		fi
	fi
	
	for i in ${FORMATS} ; do
		FORMAT="$i"
	done
	if test "${FORMAT}" != "" ; then
		GREP_FORMAT="`echo ${FORMAT} | ${SED} s/\\\\+/\\\\\\\\+/g`"
		FORMAT_INDEX=`${GREP} -n "${GREP_FORMAT}" ${BINFORMAT_FILE} | ${TAIL} -1 | ${CUT} -d: -f1`
	else
		FORMAT_INDEX=-1
		FORMAT="(unknown)"
	fi
fi
if test "$1" != "-i" ; then
	echo ${FORMAT}
else
	echo ${FORMAT_INDEX}
fi
