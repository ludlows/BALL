/*
LIGSITECS, a pocket identification method using Connolly surface
http://scoppi.biotec.tu-dresden.de

written by
   Bingding Huang
   Bioinformatics group
   Biotec & Department of Computing
   Tazberg 47, 01307
   TU Dresden, Germany
   Email bhuang@biotec.tu-dresden.de

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

#include "mySurface.h"
#include <math.h>
#include <stdlib.h>

Vector3 getMassCenter( AtomContainer&  atom_container )
{
  Vector3 mass_center = Vector3(0., 0., 0.);
  int atom_num = 0;
  
  for( AtomIterator atom_it = atom_container.beginAtom(); +atom_it; ++atom_it )
    {
      atom_num ++;
      mass_center += atom_it->getPosition();
    }
  mass_center.x = (float)mass_center.x/atom_num;
  mass_center.y = (float)mass_center.y/atom_num;    
  mass_center.z = (float)mass_center.z/atom_num;

  return mass_center;
}

Vector3 doPreTranslation(AtomContainer&  atom_container )
{
  Vector3 mc = getMassCenter( atom_container );
	TranslationProcessor tp( -mc );
	atom_container.apply( tp );
	return mc;
  
}




Vector3 getExtent( AtomContainer& atom_container )
{

  Vector3 pos;
  float maxX=0, maxY=0,maxZ=0;
  for( AtomIterator atom_it = atom_container.beginAtom(); +atom_it; ++atom_it )
    {
      pos= atom_it->getPosition();
      if (pos.x<0)
	{
	  pos.x =-pos.x;
	}
      if (pos.y<0)
	{
	  pos.y=-pos.y;
	}
      if (pos.z<0)
	{
	  pos.z=-pos.z;
	}

      if( pos.x>maxX)
	{
	  maxX=pos.x;
	}
      if (pos.y>maxY)
	{
	  maxY=pos.y;
	}
      if(pos.z>maxZ)
	{
	  maxZ=pos.z;
	}

    }
  return Vector3(maxX,maxY,maxZ);
}




  
Vector3 getGrid_origin(int GRID_SIZE, float GRID_SPACING)
{
  float r;
  int   r_idx;
  
  if(GRID_SIZE % 2 == 0)
    {
      r     = GRID_SIZE * GRID_SPACING / 2.0; // Angstrom
      r_idx = GRID_SIZE / 2; // index
    }
  else
    {
      r     = (GRID_SIZE - 1.0) * GRID_SPACING / 2.0; // Angstrom
      r_idx = (GRID_SIZE - 1) / 2; // index
    }
  
  return Vector3(r, r, r); // in unit of Angstroms
  
}


Vector3 getGrid_origin(Vector3 GRID, float GRID_SPACING)
{
  float rx,ry,rz;
  int   r_idx;
  
  if(int(GRID.x)% 2 == 0)
    {
      rx     = GRID.x * GRID_SPACING / 2.0; // Angstrom
    }
  else
    {
      rx     = (GRID.x - 1.0) * GRID_SPACING / 2.0; // Angstrom
    }
		
  if(int(GRID.y)% 2 == 0)
    {
      ry     = GRID.y * GRID_SPACING / 2.0; // Angstrom
    }
  else
    {
      ry     = (GRID.y - 1.0) * GRID_SPACING / 2.0; // Angstrom
    }
  
  if(int(GRID.z)% 2 == 0)
    {
      rz     = GRID.z * GRID_SPACING / 2.0; // Angstrom
    }
  else
    {
      rz     = (GRID.z - 1.0) * GRID_SPACING / 2.0; // Angstrom
    }
  
  return Vector3(rx, ry, rz); // in unit of Angstroms
  
}
