/*
http://projects.biotec.tu-dresden.de/pocket/download.html


"pocket.pdb", the identified pocket sites in pdb format
"pocket_all.pdb", all the pocket sites in pdb format 
"pocket_r.pdb", all the grid points before clustering in pdb format

and one python script file called "pocket.py" to visualize the pocket sites using PyMol. 
If you have pymol on your machine, just type:
pymol pocket.py 




LIGSITECS, a pocket identification method using Connolly surface
http://scoppi.biotec.tu-dresden.de

written by
   Bingding Huang
   Bioinformatics group
   Biotec & Department of Computing
   Tazberg 47, 01307
   TU Dresden, Germany
   Email bhuang@biotec.tu-dresden.de


LICENSE:
LIGSITECS, a pocket identification method using Connolly surface
Copyright (C) LIGSITECS  Bingding Huang
This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
*/

#include <iostream>
#include <set>
#include "mySurface.h"
#include <BALL/STRUCTURE/analyticalSES.h>
#include <BALL/KERNEL/system.h>
#include <BALL/FORMAT/PDBFile.h>
#include <BALL/STRUCTURE/defaultProcessors.h>
#include <BALL/DATATYPE/string.h>
#include <BALL/STRUCTURE/geometricProperties.h>
#include <BALL/STRUCTURE/geometricTransformations.h>
#include <BALL/STRUCTURE/surfaceProcessor.h>
#include <BALL/MATHS/surface.h>
#include <BALL/MATHS/vector4.h>
#include <BALL/SYSTEM/timer.h>
#include <BALL/COMMON/exception.h>
#include <BALL/STRUCTURE/numericalSAS.h>
#include <BALL/DATATYPE/hashMap.h>

#include <stdlib.h>

using namespace BALL;
using namespace std;

std::vector<std::vector<float>> getPockets(string radius_file_path, string pdb_file_name, float grid_space=1.0, int SSSthreshold=6,float density=0.5, int numberofpockets=3);
