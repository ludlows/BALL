// -*- Mode: C++; tab-width: 2: -*-
// vi: set ts=2:
//
//

#ifndef BALL_QSAR_SURFACEBASE_H
#define BALL_QSAR_SURFACEBASE_H

#ifndef BALL_CONCEPT_MOLECULE_H
#include <BALL/KERNEL/molecule.h>
#endif

#ifndef BALL_QSAR_DESCRIPTOR_H
#include <BALL/QSAR/descriptor.h>
#endif

namespace BALL
{
	/** Generic QSAR molecular surface descriptors class
		\\
	*/
	class BALL_EXPORT SurfaceBase
		:	public Descriptor
	{
		public:

		BALL_CREATE(SurfaceBase)

		/** @name Constructors and Destructors
		*/
		//@{
		/** Default constructor
		*/
		SurfaceBase();

		/** Copy constructor
		*/
		SurfaceBase(const SurfaceBase& sb);

		/** Named constructor
		*/
		SurfaceBase(const String& name);

		/** Named unity constructor
		*/
		SurfaceBase(const String& name, const String& unit);

		/** Destructor
		*/
		virtual ~SurfaceBase();
		//@}

		/** @name Assignment
		*/
		//@{
		/** Assignment operator
		*/
		virtual SurfaceBase& operator = (const SurfaceBase& sb);
		//@}

		protected:

		/** @name Predicates
		*/
		//@{
		bool isValid(Molecule& molecule);
		//@}
		
		/** @name Accessors
		*/
		//@{
		void calculate(Molecule& molecule);
		//@}
	};
} // namespace BALL

#endif

