// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// $Id: pyChainList.h,v 1.11.4.2 2005/08/19 11:25:55 amoll Exp $
//

#ifndef BALL_PYTHON_PYCHAINLIST_H
#define BALL_PYTHON_PYCHAINLIST_H

#ifndef BALL_COMMON_H
#	include <BALL/common.h>
#endif

#ifndef BALL_DATAYPE_LIST_H
#	include <BALL/DATATYPE/list.h>
#endif

namespace BALL 
{
	class Chain;
	class AtomContainer;
	
	/** Equivalent for a STL::List of Chain Pointers in Python
			\ingroup PythonExtensions
	*/
	class PyChainList
		:	public List<Chain*>
	{
		public:

		BALL_CREATE(PyChainList)

		/**	@name	Type Definitions
		*/
		//@{

		/**	Chain* type
		*/
		typedef Chain* ValueType;

		/**	Pointer type
		*/
		typedef Chain** PointerType;

		/**	Iterator type.
		*/
		typedef List<Chain*>::iterator Iterator;

		/**	Constant iterator type.
		*/
		typedef List<Chain*>::const_iterator ConstIterator;

		//@}
		/**	@name	Constructors and Destructors 
		*/
		//@{

		/**	Default constructor.
				Create an empty list.
		*/
		PyChainList();

		/** Copy constructor.
				Create a copy of an existing list.
				@param	map the list to be copied
				@param	deep ignored
		*/
		PyChainList(const PyChainList& new_list);
			
		/**	Construct from a AtomContainer
				This constructor creates an PyChainList object from
				all atoms of a  \link AtomContainer AtomContainer \endlink  object.
		*/
		PyChainList(const AtomContainer& fragment, bool selected_only = false);

		/**	Destructor
		*/
		virtual ~PyChainList() throw();

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

#endif // BALL_PYTHON_PYCHAINLIST_H
