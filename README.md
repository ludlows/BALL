# BALL 1.1.1 and Ligsite csc method for pocket position prediction
BALL project 1.1.1 bugs-free version for automatic identification of pockets on protein surface using Connolly surface and degree of conservation

the details about this algorithm is expressed in the paper: https://bmcstructbiol.biomedcentral.com/articles/10.1186/1472-6807-6-19


## Welcome to BALL-1.1.1 

This is the BALL-1.1.1 bugs-free version.

I fixed some bugs. All rights are reserved by Ball-Project team.

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

# Ligsite CSC : predicting ligand binding sites using the Connolly surface and degree of conservation
if you want to use ligsite algorithm, please visit the following link:

http://projects.biotec.tu-dresden.de/pocket/download.html 

Here is my procedure to compiling ligsite CSC code.

$source ~/.bashrc  #for linux
$source ~/.bash_profile # for Mac


$cd BALL/lcs/
$ make lcs
for mac OS, remove the dependency of -lnsl in Makefile, before run $make lcs .

# Note: this BALL code is useful only for ligsite CSC
use the latest version of BALL project, if you have interests in developing new algorithm based on BALL.  
Further information is available in BALL-Team website:
  http://www.ball-project.org
