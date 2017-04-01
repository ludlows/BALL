#include "ligsite.h"

using namspace BALL;
using namspace std;


std::vector<std::vector<float>> getPockets(string radius_file_path, string pdb_file_name, float grid_space=1.0, int SSSthreshold=6,float density=0.5, int numberofpockets=3){
    Timer timer;
    timer.start();
    PDBFile PDB (pdb_file_name);
    System sys;
    PDB >> sys;
    PDB.close();
    Vector3 massCenter=doPreTranslation(sys);

  
    Vector3 extent=getExtent(sys);
    float x=extent.x*2;
    float y=extent.y*2;
    float z=extent.z*2 ;
  
    cout<<"Grid space: "<<grid_space<<endl;
    cout<<"Number of pockets: "<<numberofpockets<<endl;
    cout<<"Density: "<<density<<endl;
    cout<<"SSS threshold: "<<SSSthreshold<<endl;

    AssignRadiusProcessor ar(radius_file_path);
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


  //  string pocketfilename="pocket_r.pdb";
  // FILE*	pocketfile= fopen( pocketfilename.c_str() , "w" ) ;

  
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
		
		
		// fprintf(pocketfile, "ATOM  %5d  Fe  PKT Z%4d    %8.3f%8.3f%8.3f\n",  index_j, it->n, pos.x, pos.y, pos.z);

      pocket_vec.push_back(p);

    }
  // fclose(pocketfile);

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
	

  index_j=0;

  ans = std::vector<std::vector<float>>

  //int numberofpockets=5;
  for( it = pocket_nr.begin(); it != pocket_nr.end(); it++ ){
  	Pocket p=*it;
    index_j++;
      
    if(index_j<=numberofpockets){
		std::vector<float> pocket_pos;
		pocket_pos.push_back((it->position).x);
		pocket_pos.push_back((it->position).y);
		pocket_pos.push_back((it->position).z);
		ans.push_back(pocket_pos);
	}
	else{
		break;
		}
    }
  timer.stop();
  float USETIME =(float) timer.getCPUTime();
  cout<<"Done!"<<endl;
  cout<<"Total used time: "<<USETIME<<" seconds"<<endl; 

  return ans;
  
 
}
int main (){
	
	string radius_file_path, pdb_file_name;
	radius_file_path = "./PARSE.siz"
	pdb_file_name = "./1dwd.pdb"
	std::vector<std::vector<float>> ans;
	ans = getPockets(radius_file_path, pdb_file_name);
	int len = ans.size();
	for (int i=0;i<len;i++){
		printf("point1 %d\n", i);
		int j = ans[i].size();
		for(int k=0;k<j;k++){
			printf("%.5f ", ans[i][j]);
		}
		printf('\n');
	}
	return 0;
}

