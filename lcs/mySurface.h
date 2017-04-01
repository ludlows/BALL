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

#include <BALL/KERNEL/system.h>

#include <BALL/STRUCTURE/numericalSAS.h>
#include <BALL/STRUCTURE/analyticalSES.h>
#include <BALL/DATATYPE/hashMap.h>
#include <BALL/DATATYPE/string.h>
#include <BALL/MATHS/vector4.h> 
#include <BALL/KERNEL/residue.h>
#include <BALL/STRUCTURE/geometricProperties.h>
#include <BALL/STRUCTURE/geometricTransformations.h>


using namespace BALL;
using namespace std;


Vector3 getExtent( AtomContainer& atom_container );
int initGridSizes( AtomContainer& atom_container, float GRID_SPACING);
Vector3 getGrid_origin(int GRID_SIZE, float GRID_SPACING);
Vector3 getGrid_origin(Vector3 GRID, float GRID_SPACING);

Vector3 getMassCenter( AtomContainer& atom_container );
Vector3 doPreTranslation(AtomContainer&  atom_container );

class grid
{

 public:
  Vector3 index;
  Vector3 position;
  int grid_index;
  bool inside;

  bool interactWithProtein;
  bool onSurface;
  
  int num_vertice;
  
 public:
  
  grid(Vector3 pos, Vector3 ind)
    {
      position=pos;
      index=ind;

      grid_index=0;
      interactWithProtein=false;

      num_vertice=0;
      onSurface=false;
    }

  grid(Vector3 pos, Vector3 ind, int in)
    {
      position=pos;
      grid_index=in;
      index=ind;

      num_vertice=0;
      interactWithProtein=false;
      onSurface=false;
    }

  ~grid()
    {

    }
  void setPosition(Vector3 pos)
  {
    position=pos;
  }

  Vector3 getPosition()
  {
    return position;
  }

  void setIndex(Vector3 pos)
  {
    index=pos;
  }

  Vector3 getIndex()
  {
    return index;
  }

  void setGrid_index(int in)
  {
    grid_index=in;

  }
  
  int getGrid_index()
  {

    return grid_index;
  }

  bool isInteractWithProtein ()
  {
    return interactWithProtein;
  }
  
  void setInteractWithProtein()
  {
    interactWithProtein=true;
  }

  bool isOnSurface()
  {
    return onSurface;
  }

  void setOnSurface()
  {
    onSurface=true;
  }
  
  void addOneVertex()
  {
    num_vertice +=1;
  }
  int getNumOfVertice()
  {
    return num_vertice;
  }

  void setNumOfVertice(int num)
  {
    num_vertice=num;
  }


};

 

class myVertex
{
 private:
  Vector3 position;
  int atom_index;
  int index;
  HashSet<int> children;
  HashSet<int> offspring;
 public:
  myVertex()
    {
      position=Vector3(0.0,0.0,0.0);
      atom_index=0;
      index=0;
      //children=new HashSet<int>;
		
    }
  myVertex(Vector3 pos,int index)
    {
      position=pos;
      index=index;
      children.insert(index);
      offspring.insert(index);
      //children=new HashSet<int>;
      //offspring=new HashSet<int>;
      
    }
  ~myVertex()
    {
      //	delete children;
    }
  
  void addChild(int c)
  {
    children.insert(c);
    offspring.insert(c);
  }
  
  void addOffspring(int o)
  {
    offspring.insert(o);
  }
  int getNumberofOffspring()
  {
    return offspring.size();
  }
  int getNumberofChildren()
  {
    return children.size();
  }
  int getIndex()
  {
    return index;
  }
  
  Vector3 getPosition()
  {
    return position;
  }
  
  HashSet<int> getChildren()
    {
      return children;
    }
 
  HashSet<int> getOffspring()
    {
      HashSet<int> tmp=offspring;
      return tmp;
    }
};


class Pocket
{
    
 public:
  Pocket()
    {
      n = 0; 
      index=0;
      position = Vector3(0., 0., 0.); 
      added=false;
    }
  Pocket(int number, Vector3 pos, int ind)
    {
      n = number; 
      position = pos;
      index=ind;
      //neighbour.add(ind);
      added=false;
    }
    
		
  ~Pocket()
    {}
  bool operator < (const Pocket& p) const
  {
    // Note: we implement "<" as ">" on purpose
    // s.t. the largest peak value is on the top of multiset
    return n > p.n;  // based on pocket number only
  }
    
  int getNumber()
  {
    return n;
  }
    
  Vector3 getPosition()
	{
		return position;
  }
  void setPosition(Vector3 pos)
	{
		position=pos;
  }
  void setNumber(int n1)
  {
    n=n1;
  }
 
  void setIndex(int n)
  {
    index=n;

  }

  int getIndex()
  {
    return index;
  }

  float getDistance()
  {
    return position.getDistance(Vector3(0.,0.,0.0));
  }

  void addNeighbour(int n)
  {
    if (n!=index)
      {
	neighbour.insert(n);
	// setAdded();
      }
  }

  HashSet<int> getNeighbour()
    {
      HashSet<int> tmp=neighbour;

      return tmp;
    }
    
  int getSize()
  {
    return neighbour.size();
      
  }
  void setAdded()
  {
    added=true;

  }
		
  bool isAdded()
  {
    return added;
  }
 public:
  int  n;
  Vector3 position;
  int index;
  HashSet<int> neighbour; 
  bool added;
 
};
