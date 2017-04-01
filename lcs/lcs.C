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

int main (int argc, char **argv)
{
	
	
  string USAGE     = "Usage: lcs -i <PDB File name for protein>  [-s grid space]  [-n Number_of_pockets] [-t SSS_threshold] [-d Surface_density] [-o output_pdb_file_name]\n";
  
  if(argc<2)
    {
      cout<<USAGE<<endl;
      exit(1);
    }
  
  Timer timer;
  timer.start();
	string pdb_file_name="";
	float grid_space=1.0;
  int numberofpockets=3;
  int SSSthreshold=6;
  float density=0.5;
  
  string SEPARATOR = "--------------------------------------------------";
	
  // command line options
  int arg_count;
  for ( arg_count = 1; arg_count < argc; ++arg_count )
    {
      if (argv[arg_count][0] == '-')
	{
	  switch (argv[arg_count][1])
	    {
				case 'i':  pdb_file_name   = string( argv[++arg_count] ); break;
				case 's':  grid_space      = atof( argv[++arg_count] );   break;
	    case 't':  SSSthreshold    = atoi( argv[++arg_count] );   break;
	    case 'd':  density         = atof( argv[++arg_count] );   break;
	    case 'n':  numberofpockets   = atoi( argv[++arg_count] );   break;	
	    case 'h': 
	      cout << USAGE << endl;
	      exit(1);
	    default:
	      cout << "invalid argument " << argv[arg_count] << endl;
	    }
	}
      else
	{
	  cout << USAGE << endl;
	  exit(1);
	}
    }
  
  PDBFile PDB (pdb_file_name);
  System sys;
  PDB >> sys;
	PDB.close();
	Vector3 massCenter=doPreTranslation(sys);
	
	//PDBFile out(output_pdb_file_name,File::out);
	//out<<sys;
	//out.close();
  
  Vector3 extent=getExtent(sys);
  float x=extent.x*2;
  float y=extent.y*2;
  float z=extent.z*2 ;
  
  cout<<"Grid space: "<<grid_space<<endl;
  cout<<"Number of pockets: "<<numberofpockets<<endl;
  cout<<"Density: "<<density<<endl;
  cout<<"SSS threshold: "<<SSSthreshold<<endl;

  char* ball="BALL";
  char* ball_path=getenv(ball); // get the BALL folder
  string radiusFile=string(ball_path)+"/data/radii/PARSE.siz"; // the radius file

  AssignRadiusProcessor ar(radiusFile);
  sys.apply(ar); //assign the radius
	
	
  int gx=(int)ceil( (x+2.0)/grid_space ) + 1;
  int gy=(int)ceil( (y+2.0)/grid_space ) + 1;
  int gz=(int)ceil( (z+2.0)/grid_space ) + 1;
  
  Vector3 gc=Vector3(gx,gy,gz);
  
  cout<<"Grid  Size: "<<gc<<endl;
	Vector3 grid_origin=getGrid_origin(gc,grid_space);
	cout<<"Grid origin: "<< -grid_origin<<endl;
	cout<<"Mass center: "<<massCenter<<endl;
  
  SurfaceProcessor sp;
  sp.setDensity(density);
  sys.apply (sp);
  
  Surface surface = sp.getSurface ();
  
  //cout << "Number of triangle " << surface.getNumberOfTriangles() << endl;
  //cout << "Number of vertice " << surface.getNumberOfVertices() << endl;
  
  multiset<class Pocket> pocket_set;
  std::vector<grid> grids;
  
  int pocket=0;
  //  grid g;
  
  int numOfGrid=0;
  
  for(int index_x=0; index_x<gx; index_x++)
    for( int index_y=0; index_y<gy; index_y++)
      for (int index_z=0; index_z<gz; index_z++)
	{
	  float pos_x=-grid_origin.x +index_x*grid_space;
	  float pos_y=-grid_origin.y +index_y*grid_space;
	  float pos_z=-grid_origin.z +index_z*grid_space;
	  Vector3 pos=Vector3(pos_x,pos_y,pos_z);
	  grid  g=grid(Vector3(pos_x,pos_y,pos_z), Vector3((float) index_x,(float) index_y,(float) index_z), numOfGrid);
	  //cout<<g.getIndex()<<endl;
	  grids.push_back(g);
	  numOfGrid++;
	}
  

  float NEAR_RADIUS=1.6;
  float WATER_RADIUS=2.0;
  
  
  cout<<"Number of grid points: "<<numOfGrid<<endl;
  
  for ( AtomIterator it = sys.beginAtom(); +it; ++it )
    {
      Vector3 atom_position = it->getPosition();  // Angstrom
      // what is the dimension we have to check for each atom
      // in the units of GRID_SPACING (index)
      int i, j, k;
      
      int lower_loop_bound_x=(int) ceil((atom_position.x+ grid_origin.x -it->getRadius() - WATER_RADIUS) /grid_space);
      int upper_loop_bound_x = (int)floor((atom_position.x+ grid_origin.x +it->getRadius() + WATER_RADIUS) /grid_space);
      
      int lower_loop_bound_y = (int) ceil((atom_position.y+ grid_origin.y -it->getRadius() - WATER_RADIUS) /grid_space);
      int upper_loop_bound_y = (int)floor((atom_position.y +grid_origin.y +it->getRadius() + WATER_RADIUS) /grid_space);
      
      int lower_loop_bound_z = (int) ceil((atom_position.z+ grid_origin.z -it->getRadius() - WATER_RADIUS) /grid_space); 
      int upper_loop_bound_z = (int)floor((atom_position.z+ grid_origin.z +it->getRadius() + WATER_RADIUS) /grid_space);
      
      // check whether the loop points are out of the global grid bound
      if ( lower_loop_bound_x < 0 ) { lower_loop_bound_x =0;}
      if ( upper_loop_bound_x > gx) { upper_loop_bound_x = gx;}
      if ( lower_loop_bound_y < 0 ) { lower_loop_bound_y = 0;}
      if ( upper_loop_bound_y > gy ) { upper_loop_bound_y =gy;}
      if ( lower_loop_bound_z < 0 ) { lower_loop_bound_z = 0;}
      if ( upper_loop_bound_z > gz ) { upper_loop_bound_z = gz;}
      
	// iterate all the grids around the atom
	// in the units of grid points (index)
      Vector3 pos;
      
      for ( i = lower_loop_bound_x; i < upper_loop_bound_x; i++ )
	for ( j = lower_loop_bound_y; j < upper_loop_bound_y; j++ )
	  for ( k = lower_loop_bound_z; k < upper_loop_bound_z; k++ )
	    {
	      // set the position we want to check, in the units of Angstrom
	      
	      pos=Vector3 (-grid_origin.x+ i *grid_space, -grid_origin.y+j*grid_space, -grid_origin.z +k *grid_space );
	      
	      int index=i*gy*gz + j*gz +k;
	      
	      if( !grids[index].isInteractWithProtein() and pos.getDistance(atom_position) <= (it->getRadius()+NEAR_RADIUS)  )
		{
		  grids[index].setInteractWithProtein();
		  
		}
	    }
    }
  
    
  //cout<<"Finish steric checking"<<endl;
  

  Vector3 pos;
  WATER_RADIUS=1.0;
  float PROBE_RADIUS=1.0;

  
  for (int v=0;v<surface.getNumberOfVertices();v++)
    {	
      
      Vector3 atom_position = surface.getVertex(v);
      
      // what is the dimension we have to check for each atom
      // in the units of GRID_SPACING (index)
      int i, j, k;
      
      int lower_loop_bound_x=(int) ceil((atom_position.x+ grid_origin.x  - WATER_RADIUS) /grid_space);
      int upper_loop_bound_x = (int)floor((atom_position.x+ grid_origin.x  + WATER_RADIUS) /grid_space);
      
      int lower_loop_bound_y = (int) ceil((atom_position.y+ grid_origin.y  - WATER_RADIUS) /grid_space);
      int upper_loop_bound_y = (int)floor((atom_position.y +grid_origin.y  + WATER_RADIUS) /grid_space);
      
      int lower_loop_bound_z = (int) ceil((atom_position.z+ grid_origin.z  - WATER_RADIUS) /grid_space); 
      int upper_loop_bound_z = (int)floor((atom_position.z+ grid_origin.z  + WATER_RADIUS) /grid_space);
      
      // check whether the loop points are out of the global grid bound
      if ( lower_loop_bound_x < 0 ) { lower_loop_bound_x =0;}
      if ( upper_loop_bound_x > gx) { upper_loop_bound_x = gx;}
      if ( lower_loop_bound_y < 0 ) { lower_loop_bound_y = 0;}
      if ( upper_loop_bound_y > gy ) { upper_loop_bound_y =gy;}
      if ( lower_loop_bound_z < 0 ) { lower_loop_bound_z = 0;}
      if ( upper_loop_bound_z > gz ) { upper_loop_bound_z = gz;}
      
	// iterate all the grids around the atom
	// in the units of grid points (index)
      Vector3 pos;
      
      for ( i = lower_loop_bound_x; i < upper_loop_bound_x; i++ )
	for ( j = lower_loop_bound_y; j < upper_loop_bound_y; j++ )
	  for ( k = lower_loop_bound_z; k < upper_loop_bound_z; k++ )
	    {
	      // set the position we want to check, in the units of Angstrom
	      
	      pos=Vector3 (-grid_origin.x+ i *grid_space, -grid_origin.y+j*grid_space, -grid_origin.z +k *grid_space );
	      //pos.set( i *grid_space,+j*grid_space, +k *grid_space );
	      
	      int index=i*gy*gz + j*gz +k;
	      
	      if ( !grids[index].isOnSurface()  and pos.getDistance(atom_position) <= PROBE_RADIUS )
		{	
		  grids[index].setInteractWithProtein();
		  grids[index].setOnSurface();
		  
		}
	      
	    }
    }
  


  
  int ix,iy,iz;
  int index_x,index_y,index_z;
  int index;

  bool touchOneSide=false;
  bool touchBothSide=false;

  int dialog_len= (int) ceil( sqrt( float (gz*gz + gy*gy + gx*gx)) ) +1;

  int pocket_num=0;
  
  int NumberOfSSS;

  
  vector<grid>::iterator it_g;    
  
  for (it_g=grids.begin();it_g!=grids.end();it_g++)
    {
      if(!it_g->isInteractWithProtein())
	{
	  grid p= grid (*it_g);
	  ix=(int) p.getIndex().x;
	  iy=(int) p.getIndex().y;
	  iz=(int) p.getIndex().z;
	  
	  
	  NumberOfSSS=0;
	  
	  touchOneSide=false;
	  touchBothSide=false;
	 
	  //printf ("Scan along x--\n ");
	  for (index_x=0; index_x<ix;index_x++) 
	    {
	      index=iz + iy*gz + index_x*gz*gy;
	      if (grids[index].isOnSurface() or grids[index].isInteractWithProtein())
		{
		  touchOneSide=true;
		  break;
		}
	    }
	 
	  if (touchOneSide)
	    {
	      //#//print "Scan along x++ "
	      for (index_x=ix+1;index_x<gx;index_x++)
		{
		  index=iz + iy*gz + index_x*gz*gy;
		  if (grids[index].isOnSurface() or grids[index].isInteractWithProtein())
		    { 
		      touchBothSide=true;
					
		      break;
		    }
		}
	    }
	  if (touchBothSide)
	    {
	      NumberOfSSS++;
	    }
        


	  //printf ("Scan along y--\n ");
	  touchOneSide=false;
	  touchBothSide=false;
        

	  for (index_y=0; index_y<iy;index_y++) 
	    {
	      index=iz + index_y*gz + ix*gz*gy;
	      if (grids[index].isOnSurface() or grids[index].isInteractWithProtein())
		{
		  touchOneSide=true;
		  break;
		}
	    }
	 
	  if (touchOneSide)
	    {
	      //#print "Scan along y++ "
	      for (index_y=iy+1;index_y<gy;index_y++)
		{

		  index=iz + index_y*gz + ix*gz*gy;
		  if (grids[index].isOnSurface() or grids[index].isInteractWithProtein())
		    { 
		      touchBothSide=true;
		      break;
		    }
		}
	    }
	  if (touchBothSide)
	    {
	      NumberOfSSS++;
	    }
	 
	 
	  //printf ("Scan along z-- \n");
	  touchOneSide=false;
	  touchBothSide=false;
        
	 
	  for (index_z=0; index_z<iz;index_z++) 
	    {
	      index=index_z + iy*gz + ix*gz*gy;
	      if (grids[index].isOnSurface() or grids[index].isInteractWithProtein())
		{
		  touchOneSide=true;
		  break;
		}
	    }
	 
	  if (touchOneSide)
	    {
	      //#print "Scan along z++ "
	      for (index_z=iz+1;index_z<gz;index_z++)
		{
		  index=index_z + iy*gz + ix*gz*gy;
		  if (grids[index].isOnSurface() or grids[index].isInteractWithProtein())
		    { 
		      touchBothSide=true;
		      break;
		    }
		}
	    }
	 
	  if (touchBothSide)
	    {
	      NumberOfSSS++;
	
	    }


        

	  int tag;
                
	  int x=ix;
	  int  y=iy;
	  int z=iz;
        
	 
	  touchOneSide=false;
	  touchBothSide=false;
	  //printf ("Scan dialog x-1, y-1, z+1\n ");

	  for (tag=0;tag<dialog_len;tag++)
	    {
	      x= x- 1;
	      y=y-1;
	      z=z + 1;
	      if( (x >= gx or x<0 )  or (y >=gy or y<0 ) or (z >=gz or z<0) ) 
		{
		  break;
		}
          
	      index=x*gy*gz+ y*gz +z;
	      if (grids[index].isOnSurface() or grids[index].isInteractWithProtein())
		{
		  touchOneSide=true;
		  break;
		}
	    }
        
	  if (touchOneSide)
	    {
	      x=ix;
	      y=iy;
	      z=iz;
	      // #print "Scan dialog x+1, y+1, z-1, "
	      for (tag=0;tag<dialog_len;tag++)
		{
		  x= x+ 1;
		  y= y + 1;
		  z= z - 1;
		  if( (x >= gx or x<0 )  or (y >=gy or y<0 ) or (z >=gz or z<0) ) 


		    {
		      break;
		    }
		
		  index=x*gy*gz+ y*gz +z;
		  if (grids[index].isOnSurface() or grids[index].isInteractWithProtein())
		    {
		      touchBothSide=true;
		      break;
		    }
		}
	    }
	  if (touchBothSide)
	    {
	      NumberOfSSS +=1;
	    }
        
        
	  x=ix;
	  y=iy;
	  z=iz;
        
	
	  touchOneSide=false;
	  touchBothSide=false;
	  //printf ("Scan dialog x-1, y-1, z-1 \n ");
	  for (tag=0;tag<dialog_len;tag++)
	    {
	      x= x- 1;
	      y=y-1;
	      z=z - 1;
	      if( (x >= gx or x<0 )  or (y >=gy or y<0 ) or (z >=gz or z<0) ) 

		{
		  break;
		}
          
	      index=x*gy*gz+ y*gz +z;
	      if (grids[index].isOnSurface() or grids[index].isInteractWithProtein())
		{
		  touchOneSide=true;
		  break;
		}
	    }
        
	  if (touchOneSide)
	    {
	      x=ix;
	      y=iy;
	      z=iz;
	      // #print "Scan dialog x+1, y+1, z+1, "
	      for (tag=0;tag<dialog_len;tag++)
		{
		  x= x+ 1;
		  y= y + 1;
		  z= z + 1;
		  if( (x >= gx or x<0 )  or (y >=gy or y<0 ) or (z >=gz or z<0) ) 

		    {
		      break;
		    }
		 
		  index=x*gy*gz+ y*gz +z;
		  if (grids[index].isOnSurface() or grids[index].isInteractWithProtein())
		    {
		      touchBothSide=true;
		      break;
		    }
		 
		}
	    }
	  if (touchBothSide)
	    {
	      NumberOfSSS +=1;
	    }
        
	

	 
	  x=ix;
	  y=iy;
	  z=iz;
        
	
	  touchOneSide=false;
	  touchBothSide=false;
	  //printf ("Scan dialog x+1, y-1, z+1 \n ");
	  for (tag=0;tag<dialog_len;tag++)
	    {
	      x= x + 1;
	      y=y - 1;
	      z=z + 1;
	      if( (x >= gx or x<0 )  or (y >=gy or y<0 ) or (z >=gz or z<0) ) 

		{
		  break;
		}
          
	      index=x*gy*gz+ y*gz +z;
	      if (grids[index].isOnSurface() or grids[index].isInteractWithProtein())
		{
		  touchOneSide=true;
		  break;
		}
	    }
        
	  if (touchOneSide)
	    {
	      x=ix;
	      y=iy;
	      z=iz;
	      // #print "Scan dialog x-1, y+1, z-1, "
	      for (tag=0;tag<dialog_len;tag++)
		{
		  x= x - 1;
		  y= y + 1;
		  z= z - 1;
		  if( (x >= gx or x<0 )  or (y >=gy or y<0 ) or (z >=gz or z<0) ) 

		    {
		      break;
		    }
		 
		  index=x*gy*gz+ y*gz +z;
		  if (grids[index].isOnSurface() or grids[index].isInteractWithProtein())
		    {
		      touchBothSide=true;
		      break;
		    }
		}
	    }
	  if (touchBothSide)
	    {
	      NumberOfSSS +=1;
	    }
       

	

 
	  x=ix;
	  y=iy;
	  z=iz;
        
	
	  touchOneSide=false;
	  touchBothSide=false;
	  //#print "Scan dialog x+1, y-1, z-1, "
	  for (tag=0;tag<dialog_len;tag++)
	    {
	      x= x + 1;
	      y=y - 1;
	      z=z  -1;
	      if( (x >= gx or x<0 )  or (y >=gy or y<0 ) or (z >=gz or z<0) ) 

		{
		  break;
		}
	     
	      index=x*gy*gz+ y*gz +z;
	      if (grids[index].isOnSurface() or grids[index].isInteractWithProtein())
		{
		  touchOneSide=true;
		  break;
		}
	    }
	 
	  if (touchOneSide)
	    {
	      x=ix;
	      y=iy;
	      z=iz;
	      // #print "Scan dialog x-1, y+1, z+1, "
	      for (tag=0;tag<dialog_len;tag++)
		{
		  x= x - 1;
		  y= y + 1;
		  z= z + 1;
		  if( (x >= gx or x<0 )  or (y >=gy or y<0 ) or (z >=gz or z<0) ) 
		    {
		      break;
		    }
		 
		  index=x*gy*gz+ y*gz +z;
		  if (grids[index].isOnSurface() or grids[index].isInteractWithProtein())
		    {
		      touchBothSide=true;
		      break;
		    }
		}
	    }
        
	  if (touchBothSide)
	    {
	      NumberOfSSS +=1;
	    }
	 



	  p.setNumOfVertice(NumberOfSSS);
	  //#print "Number of PSP", NumberOfSSS
	 
	  if (NumberOfSSS>=SSSthreshold)
	    //print "deep buried ", p
	    {
	      Pocket p= Pocket(NumberOfSSS,it_g->getPosition(),pocket_num);
	      pocket_num +=1;
	      pocket_set.insert(p);
	     
	    }
	 
	  NumberOfSSS=0;
	 
	}
      
    }
  





  cout<<"Number of redundant pockets "<<pocket_set.size()<<endl;

  multiset<class Pocket>::iterator it;//
  

  multiset<class Pocket> pocket_copy;
  int k=0;
  for(it = pocket_set.begin(); it != pocket_set.end(); it++ )
    {
      if(k<pocket_set.size())
	//if(k<500)
	{
	  pocket_copy.insert(*it);
	  k++;
	}
      
      else{
	break;
      }

    }
  
 
  multiset<class Pocket> pocket_nr;
	
  string pocketfilename="pocket_r.pdb";
  FILE*	pocketfile= fopen( pocketfilename.c_str() , "w" ) ;

  
  std::vector<Pocket> pocket_vec;

  int index_j=0;

  for( it = pocket_copy.begin(); it != pocket_copy.end(); it++ )
	{
		index_j++;
		
		Pocket p=*it;
		
		Vector3 pos=p.getPosition();
		pos =pos+massCenter; //back to oringal input protein
		p.setPosition(pos);
		pos=p.getPosition();
		
		
		fprintf(pocketfile, "ATOM  %5d  Fe  PKT Z%4d    %8.3f%8.3f%8.3f\n",  index_j, it->n, pos.x, pos.y, pos.z);

      pocket_vec.push_back(p);

    }
  fclose(pocketfile);

  int size=pocket_vec.size();


  //cout<<"Size of vec "<<size<<endl;
  //cout<<"Now caluate distance"<<endl;
  for(int m=0; m< size-1; ++m)
    {
      Vector3 pos1=pocket_vec[m].getPosition();

      float mindis=100;
      // cout<<" m "<< m<<endl;
      for (int n=1; n<size;++n)
	{
	  if (n!=m)
	    {
	      Vector3 pos2=pocket_vec[n].getPosition();
	      float d=pos1.getDistance(pos2);
	      if (d <=3.0*grid_space)
		{
		  //cout<<"Adding "<< n <<endl;
		  pocket_vec[m].addNeighbour(n);
		  //pocket_vec[n].setAdded();
		}
	    }
	}
    }	   

  
  for(int m=0; m< size; ++m)
    {
      int ns=pocket_vec[m].getSize();
      //  cout<<m<<" size: "<<ns<<endl;
      if (ns>0 && !pocket_vec[m].isAdded())
	{
	  for (int k=0; k<20;k++)
	    {
	      HashSet<int> tmp=pocket_vec[m].getNeighbour();
	      HashSet<int>::iterator it_t=tmp.begin();
	      
	      for(;it_t!=tmp.end();++it_t)
		{
		  int o= *it_t;
		  
		  if (o!=m)
		    {
		      int ns1=pocket_vec[o].getSize();
		      if (ns1>0)
			{
			  HashSet<int> tmp1=pocket_vec[o].getNeighbour();
			  HashSet<int>::iterator it_t1=tmp1.begin();
		 
			  pocket_vec[o].setAdded();
			  for(;it_t1!=tmp1.end();++it_t1)
			    {
			      int o1= *it_t1;
			      
			      if ( !pocket_vec[o1].isAdded() &&o1!=m )
				{
				  pocket_vec[m].addNeighbour(o1);
				  pocket_vec[o1].setAdded();
				}
			    }
			}
		    }
		}

	    }
	  
	}
    }




  for(int m=0; m< size; ++m)
    {
      int ns=pocket_vec[m].getSize();
      Pocket p =pocket_vec[m];
      
      if (ns>0 && !pocket_vec[m].isAdded())
	{
	  //cout<<m<<" size: "<<ns<<endl;
	  Vector3 mc=p.getPosition();
	  p.setNumber(ns+1);
	  
	  HashSet<int> tmp=pocket_vec[m].getNeighbour();
	  HashSet<int>::iterator it_t=tmp.begin();
	      
	  for(;it_t!=tmp.end();++it_t)
	    {
	      int o= *it_t;
	      mc =mc + pocket_vec[o].getPosition();
	      

	    }
	  p.setPosition(mc/(ns+1));
	  pocket_nr.insert(p);
	}

      if(ns==0)
	{
	  p.setNumber(1);
	  pocket_nr.insert(p);
	}
    }
	
  string pocketfilename1="pocket.pdb";
  FILE*	pocketfile1= fopen( pocketfilename1.c_str() , "w" ) ;

	
  string pocketfilename2="pocket_all.pdb";
  FILE*	pocketfile2= fopen( pocketfilename2.c_str() , "w" ) ;
  

  //cout<<"Number of non-redundant pockets "<<pocket_nr.size()<<endl;

  index_j=0;

  //int numberofpockets=5;
  for( it = pocket_nr.begin(); it != pocket_nr.end(); it++ )
    {
      Pocket p=*it;
      index_j++;
      
      if(index_j<=numberofpockets)
	{
	  fprintf(pocketfile1, "ATOM  %5d  Fe  PKT Z%4d    %8.3f%8.3f%8.3f\n",  index_j, it->n, (it->position).x, (it->position).y, (it->position).z);
	}
      fprintf(pocketfile2, "ATOM  %5d  Fe  PKT Z%4d    %8.3f%8.3f%8.3f\n",  index_j, it->n, (it->position).x, (it->position).y, (it->position).z);
      
    }
  
  fclose(pocketfile1);
	fclose(pocketfile2);
	FILE*	pythonFile= fopen( "pocket.py" , "w" ) ;
	fprintf(pythonFile,"#!/usr/bin/python\n#a pymol script to visulize the pocket sites\n");
	fprintf(pythonFile,"from pymol import cmd\n");
	fprintf(pythonFile,"cmd.load(\"%s\")\ncmd.color(\"green\")\ncmd.hide()\ncmd.show(\"surface\")\ncmd.load(\"pocket.pdb\")\n",pdb_file_name.c_str());
	
	fprintf(pythonFile,"cmd.set(\"sphere_scale\",0.8)\ncmd.show(\"sphere\",\"pocket\")\ncmd.color(\"red\",\"pocket\")\n");
	//fprintf(pythonFile,"# select the atoms involved in the first pocket site\n");
	//fprintf(pythonFile,"cmd.select(\"bs\",\"(index 1 and pocket) around 5.0 and (not pocket)\")\n");
	//fprintf(pythonFile,"cmd.color(\"blue\",\"bs\")\ncmd.save(\"bs.pdb\",\"bs\")\n)";
	fprintf(pythonFile,	"cmd.deselect()\ncmd.zoom()\n");
	
	fclose(pythonFile);
  
  timer.stop();
  float USETIME =(float) timer.getCPUTime();
  cout<<"Done!"<<endl;
  cout<<"Total used time: "<<USETIME<<" seconds"<<endl; 
  
}
