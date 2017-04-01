// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// $Id: numericalSAS.h,v 1.27.4.1 2005/07/29 12:38:09 amoll Exp $
//

#ifndef BALL_STRUCTURE_NUMERICALSAS_H
#define BALL_STRUCTURE_NUMERICALSAS_H

#ifndef BALL_COMMON_H
#	include <BALL/common.h>
#endif

#ifndef BALL_MATHS_VECTOR3_H
#	include <BALL/MATHS/vector3.h>
#endif

#ifndef BALL_MATHS_SURFACE_H
#	include <BALL/MATHS/surface.h>
#endif

namespace BALL 
{
	class Atom;
	class AtomContainer;
	template <typename Key, typename T>
  class HashMap;

	/**	@name	Fast Numerical Calculation of Solvent Accessible Surface Area.
			These functions use the algorithm by Eisenhaber, Lijnzaad, Argos, Sander,
			and Scharf ("The Double Cubic Lattice Method: Efficient Approaches to 
			numerical Integration of Surface Area and Volume and to Dot
			Surface Contouring of Molecular Assemblies", J. Comput. Chem. (1995),
			<b>  15 </b>, 273-284).
	
			\ingroup Surface
			@{
	*/

	
	/** Calculate the solvent accessible surface area numerically.
			This method returns the total 
			Solvent Accessible Surface (SAS) of a BALL kernel object. 
			Atoms with a radius of 0 are ignored.
			@param  fragment the kernel object containing the atoms
			@param  probe_radius the probe radius used for the SAS
			@param  number_of_dots the number of dots used per atom
			@return the total SAS area in $ A^2$
	*/
	BALL_EXPORT float calculateSASArea(const AtomContainer& fragment, float probe_radius = 1.5,
												 Size number_of_dots = 400); 

	/** Calculate the solvent accessible volume numerically.
			This method returns the total volume enclosd by the 
			Solvent Accessible Surface (SAS) of a BALL kernel object. 
			Atoms with a radius of 0 are ignored.
			@param  fragment the kernel object containing the atoms
			@param  probe_radius the probe radius used for the SAS
			@param  number_of_dots the number of dots used per atom
			@return the volume in $ A^3$
	*/
	BALL_EXPORT float calculateSASVolume(const AtomContainer& fragment, float probe_radius = 1.5,
													 Size number_of_dots = 400); 

	/**	Calculate the Solvent Accessible Surface area for each atom.
			This method returns the surface fraction of each atom at the 
			Solvent Accessible Surface (SAS). Atoms with a radius of 0 are 
			ignored. All areas are in $ A^2$.
			@param  atom_areas a hash map containing the areas of the atoms (returned)
			@param	fragment the kernel object containing the atoms
			@param	probe_radius the probe radius used for the SAS
			@param	number_of_dots the number of dots used per atom
			@return the total SAS area in $ A^2$
	*/
	BALL_EXPORT float calculateSASAtomAreas(const AtomContainer& fragment, HashMap<const Atom*,float>& atom_areas,
															float probe_radius = 1.5, Size number_of_dots = 400);
	
	/**	Calculate a point set on the Solvent Accessible Surface.
			This method returns the points on the Solvent Accessible
			Surface (SAS) used to calculate the surface area.
			The  \link Surface Surface \endlink  object holds just the vertices, it
			does not contain any triangles. The normals for each point
			normals to the SAS in that point, their length equals the 
			fraction of the surface area represented by this point in $ A^2$:
			\f[	
					|\vec{n_i}| = \frac{\mathrm{SAS\quad\ of \quad atom} \quad i}{\mathrm{number\quad of\quad points\quad on\quad the\quad SAS\quad of\quad atom} \quad i}     
			\f]
			Atoms with a radius of 0 are ignored.
			@param  surface_points a surface object containing the point coordinates and their normals (returned)
			@param	fragment the kernel object containing the atoms
			@param	probe_radius the probe radius used for the SAS
			@param	number_of_dots the number of dots used per atom
			@return the total SAS area in $ A^2$
	*/
	BALL_EXPORT float calculateSASPoints(const AtomContainer& fragment, Surface& surface_points,
													 float probe_radius = 1.5,  Size number_of_dots = 400);

	/** Calculate a point set on the Solvent Accessible Surface for each
			atom. This method returns the point sets on the SAS used to calculate
			the surface area for each atom. 
			@see calculateSASAtomAreas
			@param	fragment the kernel object containing the atoms
			@param  atom_surfaces a hashmap of atoms and Surface objects containing the point sets for each atom (returned)
			@param	probe_radius the probe radius used for the SAS
			@param	number_of_dots the number of dots used per atom
			@return the total SAS area in $ A^2$
	 */
	BALL_EXPORT float calculateSASAtomPoints(const AtomContainer& fragment, 
											 				 std::vector< std::pair<Vector3, Surface> >& atom_surfaces,
															 float probe_radius = 1.5,  Size number_of_dots = 400);
   /** @} */
} // namespace BALL

#endif // BALL_STRUCTURE_NUMERICALSAS_H
