// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// $Id: geometricProperties.h,v 1.24.6.2 2005/08/13 15:58:23 amoll Exp $
//

#ifndef BALL_STRUCTURE_GEOMETRICPROPERTIES_H
#define BALL_STRUCTURE_GEOMETRICPROPERTIES_H

#ifndef BALL_COMMON_H
#	include <BALL/common.h>
#endif

#ifndef BALL_MATHS_VECTOR3_H
#	include <BALL/MATHS/vector3.h>
#endif

#ifndef BALL_MATHS_SIMPLEBOX3_H
#	include <BALL/MATHS/simpleBox3.h>
#endif

#ifndef BALL_KERNEL_ATOM_H
#	include <BALL/KERNEL/atom.h>
#endif

#ifndef BALL_KERNEL_FRAGMENT_H
#	include <BALL/KERNEL/fragment.h>
#endif

#ifndef BALL_CONCEPT_PROCESSOR_H
#	include <BALL/CONCEPT/processor.h>
#endif

#ifndef BALL_DATATYPE_STRING_H
#	include <BALL/DATATYPE/string.h>
#endif

#include <vector>

namespace BALL 
{

	/**	@name	Geometric property processors.
			The applicators, processors, and collectors described in 
			this chapter are used to extract geometric properties out
			of a given molecular object or to extract parts of these
			objects according to their geometric properties. \par
			Using the  \link BoundingBoxProcessor BoundingBoxProcessor \endlink , the bounding box 
			of a given molecular object can be calculated. The bounding box is
			represented by the lowest and highest coordinates occuring in the 
			molecular object, i.e. the bounding box is the smallest rectangular
			box (with sides parallel to the coordinate axes) that encloses all
			atoms in the molecular object. \par
			The  \link GeometricCenterProcessor GeometricCenterProcessor \endlink  calculates the geometric 
			center of all atoms contained in the molecular object it is applied to. \par
			With the aid of the  \link FragmentDistanceCollector FragmentDistanceCollector \endlink  it is possible
			to collect all molecular fragments that are within a given distance
			from a certain fragment. This is useful to extract the relevant molecular
			environment (e.g. to examin a binding site). \par
	*/
	//@{

	/**	Bounding box creating processor.
			This class iterates over all atoms of a given molecular object and
			determines the lowest and the highest coordinates occuring. It returns
			two coordinates ( \link getLower getLower \endlink ,  \link getUpper getUpper \endlink ) describing the smallest
			cuboid (whose sides are parallel to the planes defined by the corrdinate
			axes) enclosing all atoms of the molecular object. \par
			This processor is useful to determine the extent of a molecular object
			if you want to define a  \link THashGrid THashGrid \endlink  or alike objects. \par
			The coordinates returned by  \link getLower getLower \endlink  and  \link getUpper getUpper \endlink  are only
			valid, if the processor has been applied to a molecular object containing
			atoms. \par
			
			 \par
	\ingroup StructureMiscellaneous
	*/
	class BALL_EXPORT BoundingBoxProcessor
		:	public UnaryProcessor<Atom>
	{
		public:

		/** @name Processor related methods.
		*/
		//@{
			
		/**
		*/
		virtual bool start()
			throw();

		/**
		*/
		virtual bool finish()
			throw();

		/**
		*/
		virtual Processor::Result operator () (Atom& atom)
			throw() { return operator() (atom.getPosition());}

		/**
		*/
		virtual Processor::Result operator () (const Vector3& v)
			throw();


		//@}
		/**	@name Accessors
		*/
		//@{

		/** Return the bounding box
		*/
		SimpleBox3 getBox() const
			throw();

		/**	Returns the lower corner of the bounding box
		*/
		const Vector3& getLower() const
			throw();

		/**	Returns the upper corner of the bounding box
		*/
		const Vector3& getUpper() const
			throw();

		//@}
			
		private:

		Vector3	lower_;
		Vector3	upper_;
	};

	/**	Calculates the geometric center of a given Composite object.
			This processor calculates the geometric center of the atom coordinates
			of a given molecular object. \par
			The geometric center is calculated as follows: \par
			\[
				\vec{C} = \frac{1}{N} \sum_{i}{N} \vec{r_i}
			\]
			Where $\vec{r_i}$ represents the coordinates of the ith atom. \par
			\ingroup StructureMiscellaneous
	*/
	class BALL_EXPORT GeometricCenterProcessor
		:	public UnaryProcessor<Atom> 
	{
		public:

		/**	@name Processor related methods
		*/
		//@{

		/**
		*/
		virtual bool start()
			throw();

		/**
		*/
		virtual bool finish()
			throw();

		/**
		*/
		virtual Processor::Result operator()(Atom& atom)
			throw() { return operator()(atom.getPosition());}

		/**
		*/
		virtual Processor::Result operator()(const Vector3& v)
			throw();

		//@}
		/**@name	Accessors
		*/
		//@{

		/**	Returns the center of the object
		*/
		Vector3& getCenter()
			throw();

		//@}

		private:

		Vector3	center_;
		Size		n_;
	};


	/**	Collects all MolecularFragments that are close enough to another 
			molecular fragment.
			This processor examines the distances between every atom of a given fragment
			(further referred to as the reference fragment) and all other atoms in a molecular 
			object he is applied to. If any atom of a fragment is closer to any atom of the
			reference fragment, the whole fragment is collected in an array. \par
			The reference fragment itself is also contained in this array, if it is part
			of the molecular object the collector is applied to. \par
			The array only contains pointers to the fragments, the fragments are neither 
			changed, nor removed from the molecular object. \par
			The reference fragment may either be given by a specialized constructor (also
			together with the distance) or using  \link setFragment setFragment \endlink . \par
			The fragment array is emptied prior to each collection run. \par
			
			 \par
	\ingroup StructureMiscellaneous
	*/
	class BALL_EXPORT FragmentDistanceCollector
		: public UnaryProcessor<Composite> 
	{		
		public:

		/** @name Constructors and Destructors 
		*/
		//@{

		/**	Default constructor
		*/
		FragmentDistanceCollector()
			throw();

		/**	Constructor.
				Creates a new collector and sets the reference composite
				@param	composite the reference composite
		*/
		FragmentDistanceCollector(const Composite& composite)
			throw();

		/**	Constructor.
				Creates a new collector and sets the reference composite and the distance.
				@param	composite the reference composite
				@param	distance the maximum distance between any two atoms
		*/
		FragmentDistanceCollector(const Composite& composite, float distance)
			throw();
			
		virtual ~FragmentDistanceCollector()
			throw()
		{}	
			
		//@}
		/**	@name	Processor related methods
		*/
		//@{

		/**
		*/
		virtual bool start()
			throw();

		/**
		*/
		virtual bool finish()
			throw();

		/**
		*/
		virtual Processor::Result operator()(Composite& composite)
			throw();

		//@}
		/**	@name Accessors
		*/
		//@{

		/**	Returns the number of molecular fragments found
				@return	the number of fragments in the array
		*/
		Size getNumberOfFragments()
			throw();

		/**	Sets the reference composite
				@param	composite the new reference composite
		*/
		void setComposite(const Composite& composite)
			throw();

		/**	Gets the reference composite
				@return a const pointer to the reference composite
		*/
		const Composite* getComposite() const
			throw();

		/**	Gets the maximum distance
				@return the maximum distance
		*/
		float getDistance() const
			throw();
		
		/**	Sets the maximum distance
				@param	distance the new maximum distance 
		*/
		void setDistance(float distance)
			throw();

		//@}
		
		/**	The array containing all molecular fragments collected
		*/
		std::vector<Fragment*>	fragments;


		protected:

		std::vector<Fragment*>	all_fragments_;
		const Composite*	reference_composite_;
		float							squared_distance_;
	};


	//@}
	/**	@name	Angle Calculation
	\ingroup StructureMiscellaneous
	*/
	//@{
		
	/**	Calculate the torsion angle between four atoms
	*/
	BALL_EXPORT Angle calculateTorsionAngle(const Atom& a1, const Atom& a2, const Atom& a3, const Atom& a4)
		throw(Exception::IllegalPosition);

	/**	Calculate the bond angle between three atoms
	*/
	BALL_EXPORT Angle calculateBondAngle(const Atom& a1, const Atom& a2, const Atom& a3)
		throw(Exception::IllegalPosition);

	//@}
} // namespace BALL

#endif // BALL_STRUCTURE_GEOMETRICPROPERTIES_H
