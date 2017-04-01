// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// $Id: pyResidueList.h,v 1.12.4.2 2005/08/19 11:25:57 amoll Exp $
//

#ifndef BALL_PYTHON_PYRESIDUELIST_H
#define BALL_PYTHON_PYRESIDUELIST_H

#ifndef BALL_COMMON_H
#	include <BALL/common.h>
#endif

#ifndef BALL_DATAYPE_LIST_H
#	include <BALL/DATATYPE/list.h>
#endif

namespace BALL 
{
	class Residue;
	class AtomContainer;
	
	/** Equivalent for a STL::List of Residue Pointers in Python
			\ingroup PythonExtensions
	*/
	class PyResidueList
		:	public List<Residue*>
	{
		public:

		BALL_CREATE(PyResidueList)

		/**	@name	Type Definitions
		*/
		//@{

		/**	Residue* type
		*/
		typedef Residue* ValueType;

		/**	Pointer type
		*/
		typedef Residue** PointerType;

		/**	Iterator type.
		*/
		typedef List<Residue*>::iterator Iterator;

		/**	Constant iterator type.
		*/
		typedef List<Residue*>::const_iterator ConstIterator;

		//@}

		/**	@name	Constructors and Destructors */
		//@{

		/**	Default constructor.
				Create an empty list.
		*/
		PyResidueList();

		/** Copy constructor.
				Create a copy of an existing list.
				@param	map the list to be copied
				@param	deep ignored
		*/
		PyResidueList(const PyResidueList& new_list);
			
		/**	Construct from a AtomContainer
				This constructor creates an PyResidueList object from
				all residues of a  \link AtomContainer AtomContainer \endlink  object.
		*/
		PyResidueList(const AtomContainer& fragment, bool selected_only = false);

		/**	Destructor
		*/
		virtual ~PyResidueList() throw();
		//@}

		/**	@name Assignment
		*/
		//@{
		/**
		*/
		void set(const AtomContainer& fragment, bool selected_only = false);
		//@}
	};
   
} // namespace BALL

#endif // BALL_PYTHON_PYRESIDUELIST_H
