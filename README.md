# BALL
BALL project 1.1.1 bugs-free version

## Welcome to BALL-1.1.1

This is the BALL-1.1.1 bugs-free version.

I fixed some bugs. All rights are reserved by Ball-Project teem.

## Install 

$cd BALL/source

run confiure to see which software is not installed.

$./configure --disable-VIEW 

$make

$make install

edit your path as the following:

export BALL = your\_ball\_folder\_path
export LD\_LIBRARY\_PATH = \$BALL/lib/Linux\_g++\_(your version)
export PATH = "\$BALL:\$LD\_LIBRARY\_PATH:\$PATH"


if you want to use ligsite algorithm, please used the following link.

http://projects.biotec.tu-dresden.de/pocket/download.html for ligsite algorithm


Further information is available from team website:
  http://www.ball-project.org