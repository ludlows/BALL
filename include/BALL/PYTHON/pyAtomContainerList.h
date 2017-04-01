// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// $Id: pyAtomContainerList.h,v 1.10.4.2 2005/08/19 11:25:54 amoll Exp $
//

#ifndef BALL_PYTHON_PYATOMCONTAINERLIST_H
#define BALL_PYTHON_PYATOMCONTAINERLIST_H

#ifndef BALL_COMMON_H
#	include <BALL/common.h>
#endif

#ifndef BALL_DATAYPE_LIST_H
#	include <BALL/DATATYPE/list.h>
#endif

namespace BALL 
{
	class AtomContainer;
	
	/** Equivalent for a STL::List of AtomContainer Pointers in Python
			\ingroup PythonExtensions
	*/
	class PyAtomContainerList
		:	public List<AtomContainer*>
	{
		public:

		BALL_CREATE(PyAtomContainerList)

		/**	@name	Type Definitions
		*/
		//@{

		/**	AtomContainer* type
		*/
		typedef AtomContainer* ValueType;

		/**	Pointer type
		*/
		typedef AtomContainer** PointerType;

		/**	Iterator type.
		*/
		typedef List<AtomContainer*>::iterator Iterator;

		/**	Constant iterator type.
		*/
		typedef List<AtomContainer*>::const_iterator ConstIterator;

		//@}
		/**	@name	Constructors and Destructors 
		*/
		//@{

		/**	Default constructor.
				Create an empty list.
		*/
		PyAtomContainerList();

		/** Copy constructor.
				Create a copy of an existing list.
				@param	map the list to be copied
				@param	deep ignored
		*/
		PyAtomContainerList(const PyAtomContainerList& new_list);
			
		/**	Construct from a AtomContainer
				This constructor creates an PyAtomContainerList object from
				all base fragments of a  \link AtomContainer AtomContainer \endlink  object.
		*/
		PyAtomContainerList(const AtomContainer& fragment, bool selected_only = false);

		/**	Destructor
		*/
		virtual ~PyAtomContainerList() throw();
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

#endif // BALL_PYTHON_PYATOMCONTAINERLIST_H
