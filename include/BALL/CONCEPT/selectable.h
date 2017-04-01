// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// $Id: selectable.h,v 1.19.6.2 2005/07/28 13:50:52 amoll Exp $
//

#ifndef BALL_CONCEPT_SELECTABLE_H
#define BALL_CONCEPT_SELECTABLE_H

#ifndef BALL_COMMON_H
#	include <BALL/common.h>
#endif

#ifndef BALL_CONCEPT_PERSISTENTOBJECT_H
#	include <BALL/CONCEPT/persistentObject.h>
#endif

#define BALL_SELECTABLE_DEFAULT_SELECTION   false

namespace BALL 
{

	/** Selectable Concept.
			Selectable implements the ability of objects to be selected.
			It is used in the BALL kernel to mark parts of the kernel
			objects for special operations. In principle, Selectable
			provides a boolean flag and defines an interface to access
			this flag. \par

			<b>Interface:</b> Storable
			 \par
			
			 \par
			 \ingroup ConceptsMiscellaneous
	*/
	class BALL_EXPORT Selectable
	{
		public:

		/**	@name	Constructors and Destructors 
		*/
		//@{

		/**	Default Constructor.
				Creates a new selectable object and sets its state to unselected.
		*/
		Selectable()
			throw();
			
		/**	Copy constructor.
				Creates a copy of selectable object.
				@param	selectable the Selectable object to be copied
				@param	deep ignored
		*/
		Selectable(const Selectable& selectable, bool deep = true)
			throw();
			
		/**	Destructor.
				The destructor has no functionality.
		*/
		virtual ~Selectable()
			throw();

		/**	Clear the selection flag.
				Clear resets the selection flag to unselected.
		*/
		virtual void clear()
			throw();
		
		//@}

		/**	@name	Assignment 
		*/
		//@{

		/**	Assign the contents of another selectable object to this object.
				@param selectable the object whose contents are to be copied
				@param deep ignored
		*/
		void set(const Selectable& selectable, bool deep = true)
			throw();

		/**	Assignment operator.
				Assigns the contents of another Selectable object
				to this object.
				@param selectable the object to be copied
		*/
		const Selectable& operator = (const Selectable& selectable)
			throw();
		
		/**	Copy the contents of this object into another.
				@param selectable the object to be assigned to
				@param deep ignored
		*/
		void get(Selectable& selectable, bool deep = true) const
			throw();

		/**	Swap the contents of two objects.
				@param selectable the object to swap contents with
		*/
		void swap(Selectable& selectable)
			throw();
		//@}

		/**	Accessors
		*/
		//@{
		
		/**	Select the object.
				The internal flag is set to <b>true</b>.
		*/
		virtual void select()
			throw();

		/**	Deselect the object.
				The internal flag is set to <b>false</b>.
		*/
		virtual void deselect()
			throw();

		/** For faster access 
		 */
		virtual void setSelected(bool selected)
			throw() { selected_ = selected; }

		//@}
		/**	@name Predicates 
		*/
		//@{

		/**	Get the object state.
				@return bool <b>true</b>, if the object is selected, <b>false</b> otherwise
		*/
		bool isSelected() const
			throw();

		/** Equality operator
		*/
		bool operator == (const Selectable& selectable) const 
			throw();

		/** Inequality operator
		*/
		bool operator != (const Selectable& selectable) const 
			throw();

		//@}

		/**	@name Storable interface.
		*/
		//@{

		/**	Persistent stream writing.
				This method writes a boolean variable to the
				persistent stream using the <tt>writePrimitive</tt> method
				of the PersistenceManager.
				@param pm the persistence manager
		*/
		void write(PersistenceManager& pm) const
			throw();

		/** Persistent stream reading.
				This method reads a boolean variable from the
				persistent stream using the <tt>readPrimitive</tt> method
				of the PersistenceManager.
				@param pm the persistence manager
		*/
		bool read(PersistenceManager& pm)
			throw();

		//@}
		/**	@name Debugging and Diagnostics 
		*/
		//@{
		
		///
		virtual void dump(::std::ostream& s = std::cout, Size depth = 0) const
			throw();

		//@}
	
		protected:

		bool selected_;
	};

#	ifndef BALL_NO_INLINE_FUNCTIONS
#		include <BALL/CONCEPT/selectable.iC>
#	endif

} // namespace BALL

#endif // BALL_CONCEPT_SELECTABLE_H
