// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// $Id: molecule.h,v 1.34.4.2 2005/10/11 11:51:00 oliver Exp $
//

#ifndef BALL_KERNEL_MOLECULE_H
#define BALL_KERNEL_MOLECULE_H

#ifndef BALL_KERNEL_ATOMCONTAINER_H
#	include <BALL/KERNEL/atomContainer.h>
#endif

#ifndef BALL_KERNEL_FRAGMENT_H
#	include <BALL/KERNEL/fragment.h>
#endif

#ifndef BALL_KERNEL_FRAGMENTITERATOR_H
#	include <BALL/KERNEL/fragmentIterator.h>
#endif

namespace BALL 
{
	class System;

	/**	Molecule class.
			Used to represent general molecules without specific properties. \par
			
   		\ingroup KernelContainers 
	*/
	class BALL_EXPORT Molecule
		: public AtomContainer
	{
		public:

		BALL_CREATE_DEEP(Molecule)

		/**	@name	Enums
		*/
		//@{

		/**	Properties
		*/
		enum Property
		{
			IS_SOLVENT = AtomContainer::NUMBER_OF_PROPERTIES,
			NUMBER_OF_PROPERTIES
		};

		//@}
		/**	@name	Constructors and Destructors 
		*/
		//@{

		/**	Default constructor. 
		*/
		Molecule()
			throw();
	
		/** Copy constructor. 
		*/
		Molecule(const Molecule& molecule, bool deep = true)
			throw();
	
		/** Detailled constructor. 
		*/
		Molecule(const String& name)
			throw();

		/** Destructor. 
		*/
		virtual ~Molecule()
			throw();
		
		//@}
		/** @name Persistence 
		*/
		//@{

		/**	Writes a Molecule object to a persistent stream.
				@param pm the persistence manager
		*/
		void persistentWrite(PersistenceManager& pm, const char* name = 0) const
			throw(Exception::GeneralException);

		/**	Reads a Molecule object from a persistent stream.
				@param pm the persistence manager
		*/
		void persistentRead(PersistenceManager& pm)
			throw(Exception::GeneralException);

		//@}
		/**	@name	Assignment 
		*/
		//@{

		/**	Assign from another Molecule.
				@param molecule	the Molecule object to assign from
				@param  deep make a deep (=<tt>true</tt>) or shallow (=<tt>false</tt>) copy
		*/
		void set(const Molecule& molecule, bool deep = true)
			throw();

		/**	Assignment operator.
				@param molecule the Molecule to assign from
		**/
		Molecule& operator = (const Molecule& molecule)
			throw();

		/**	Assign to another Molecule.
				@param molecule	the Molecule to be assigned to
				@param  deep make a deep (=<tt>true</tt>) or shallow (=<tt>false</tt>) copy
		*/
		void get(Molecule& molecule, bool deep = true) const
			throw();

		/**	Swap the contents of two molecules.
				@param	molecule the Molecule to swap contents with
		*/
		void swap(Molecule& molecule)
			throw();
	
		//@}
		/**	@name	Accessors 
		*/
		//@{

		/**	Access the parent System.
				@return	System* pointer to the parent System
		*/
		System* getSystem()
			throw();

		/**	Get a const pointer to the parent System.
				@return	System* pointer to the parent System
		*/
		const System* getSystem() const
			throw();

		/** Insert an atom as the first child.
				@param atom the atom to add
		*/
		void prepend(Atom& atom)
			throw();

		/** Insert an atom as the last child.
				@param atom the atom to add
		*/
		void append(Atom& atom)
			throw();

		/** Insert an atom as the last child.
				@param atom the atom to add
		*/
		void insert(Atom& atom)
			throw();

		/** Insert an atom before a Composite object.
				@param atom the atom to insert
				@param before the Composite object to insert before
		*/
		void insertBefore(Atom& atom, Composite& before)
			throw();

		/** Insert an atom after a Composite object.
				@param atom the atom to insert
				@param after the Composite object to insert after
		*/
		void insertAfter(Atom& atom, Composite& after)
			throw();

		/** Remove an atom.
				@param atom the atom to remove
		*/
		bool remove(Atom& atom)
			throw();

		/** Insert an AtomContainer as the first child.
				@param atom_container the AtomContainer to add
		*/
		void prepend(AtomContainer& atom_container)
			throw();

		/** Append an AtomContainer as the last child.
				@param atom_container the AtomContainer to add
		*/
		void append(AtomContainer& atom_container)
			throw();

		/** Insert an AtomContainer as the last child.
				@param atom_container the AtomContainer to add
		*/
		void insert(AtomContainer& atom_container)
			throw();

		/** Insert an AtomContainer before a given Composite object.
				@param atom_container the AtomContainer to insert
				@param before the Composite object to insert before
		*/
		void insertBefore(AtomContainer& atom_container, Composite& before)
			throw();

		/** Insert an AtomContainer after a given Composite object.
				@param atom_container the AtomContainer to insert
				@param after the Composite object to insert after
		*/
		void insertAfter(AtomContainer& atom_container, Composite& after)
			throw();

		/**	Cut all children of <tt>atom_container</tt> and prepend them before the children of this molecule.
				@param atom_container the AtomContainer to access
		*/
		void spliceBefore(AtomContainer& atom_container)
			throw();

		/**	Cut all children of <tt>atom_container</tt> and append them after the children of this molecule.
				@param atom_container the AtomContainer to access
		*/
		void spliceAfter(AtomContainer& atom_container)
			throw();

		/**	Move the children of atom_container into this molecule.
				The children of <tt>atom_container</tt> are inserted at the position of 
				<tt>atom_container</tt> if it is a child of <tt>this</tt>.
				Otherwise the children are inserted using  \link spliceBefore spliceBefore \endlink .
		*/
		void splice(AtomContainer& atom_container)
			throw();

		/** Remove an AtomContainer.
				@param atom_container the AtomContainer to remove
		*/
		bool remove(AtomContainer& atom_container)
			throw();

		//@}
		/**	@name Debugging and Diagnostics 
		*/
		//@{

		/** Internal state dump.
				Dump the current internal state to the 
				output ostream <b>  s </b> with dumping depth <b>  depth </b>.
				@param	s output stream where to output the internal state
				@param  depth the dumping depth
		*/
		virtual void dump(std::ostream& s = std::cout, Size depth = 0) const
			throw();

		//@}

		/**	Equality operator.
				Two molecules are equal if they have the same handle.
				@see Object::operator ==.
		*/
		bool operator == (const Molecule& molecule) const
			throw();

		/**	Inequality operator
				@see operator ==
		*/
		bool operator != (const Molecule& molecule) const
			throw();


		BALL_KERNEL_DEFINE_ITERATOR_CREATORS(Fragment)


		protected:

		Molecule* getMolecule()
			throw();

		const Molecule* getMolecule() const
			throw();

		AtomContainer* getSuperAtomContainer()
			throw();

		const AtomContainer* getSuperAtomContainer() const
			throw();

		void prepend(Molecule& molecule)
			throw();

		void append(Molecule& molecule)
			throw();

		void insert(Molecule& molecule)
			throw();

		void insertBefore(Molecule& molecule, Composite& composite)
			throw();

		void insertAfter(Molecule& molecule, Composite& composite)
			throw();

		bool remove(Molecule& molecule)
			throw();

		bool isSubAtomContainerOf(const AtomContainer& atom_container) const
			throw();
	};
} // namespace BALL 

#endif // BALL_KERNEL_MOLECULE_H
