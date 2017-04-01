// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// $Id: chain.h,v 1.38.4.1 2005/07/28 14:01:54 amoll Exp $
//

#ifndef BALL_KERNEL_CHAIN_H
#define BALL_KERNEL_CHAIN_H

#ifndef BALL_KERNEL_RESIDUE_H
#	include <BALL/KERNEL/residue.h>
#endif

#ifndef BALL_KERNEL_SECONDARYSTRUCTUREITERATOR_H
#	include <BALL/KERNEL/secondaryStructureIterator.h>
#endif


#define BALL_CHAIN_DEFAULT_NAME   ' '

namespace BALL 
{
	class Protein;

	/** Protein chain class.
			This class represents a polypeptide chain within a  \link Protein Protein \endlink .
			Chains can contain  \link SecondaryStructure SecondaryStructure \endlink s or  \link Residue Residue \endlink s.
			 \par
			
	 		\ingroup KernelContainers 
	*/
	class BALL_EXPORT Chain
		: public AtomContainer
	{
		public:

		BALL_CREATE_DEEP(Chain)
		
		/**	@name	Enums
		*/
		//@{

		/**
		*/
		enum Property
		{
			NUMBER_OF_PROPERTIES = AtomContainer::NUMBER_OF_PROPERTIES
		};

		//@}
		/**	@name	Constructors and Destructors 
		*/
		//@{

		/// Default constrcutor
		Chain()
			throw();
	
		/// Copy constructor
		Chain(const Chain& chain, bool deep = true)
			throw();
	
		/// Detailled constructor
		Chain(const String& name)
			throw();

		/// Destructor
		virtual ~Chain()
			throw();
			
		//@}
		/** @name Persistence 
		*/
		//@{

		/**	Writes a Chain object to a persistent stream.
				@param pm the persistence manager
		*/
		void persistentWrite(PersistenceManager& pm, const char* name = 0) const
			throw(Exception::GeneralException);

		/**	Reads a Chain object from a persistent stream.
				@param pm the persistence manager
		*/
		void persistentRead(PersistenceManager& pm)
			throw(Exception::GeneralException);

		//@}
		/**	@name Assignment 
		*/
		//@{

		/** Assignment with cloning facility.
				The assignment is either deep or shallow (default).
				@param  chain the chain to be copied (cloned)
				@param  deep make a deep (=<tt>true</tt>) or shallow (=<tt>false</tt>) copy
		*/
		void set(const Chain& chain, bool deep = true)
			throw();

		/** Assignment operator.
				The assignment is deep.
				@param   chain the chain to be copied (cloned)
				@return  chain& - this instance.
				@see     chain::set
		*/
		Chain& operator = (const Chain& chain)
			throw();

		/** Copying with cloning facility.
				The assignment is either deep or shallow (default).
				@param  chain the chain to be assigned to
				@see    chain::set
		*/
		void get(Chain& chain, bool deep = true) const
			throw();

		/** Swapping of chains.
				@param  chain the chain this instance is being swapped with
		*/
		void swap(Chain& chain)
			throw();

		//@}
	
		/**	Equality operator.
				Two chains are equal if they have the same handle.
				@see Object::operator ==
		*/
		bool operator == (const Chain& chain) const
			throw();

		/**	Inequality operator
				@see operator ==
		*/
		bool operator != (const Chain& chain) const
			throw();


		/**	@name	Accessors 
		*/
		//@{

		/** Get a pointer to the parent protein.
				The pointer is 0 if this instance does not have a parent protein.
				@return  Protein* - mutable pointer to the parent protein
		*/
		Protein* getProtein()
			throw();

		/** Get a constant pointer to the parent protein.
				The pointer is 0 if this instance does not have a parent protein.
				@return  Protein* - constant pointer to the parent protein
		*/
		const Protein* getProtein() const
			throw();

		/** Get a pointer to a child SecondaryStructure at a given position.
				The pointer is 0 if this instance does not have a SecondaryStructure at this position.
				@param   position the position of the child SecondaryStructure
				@return  SecondaryStructure* -
								 mutable pointer to the child SecondaryStructure at <b>  position </b>
		*/
		SecondaryStructure* getSecondaryStructure(Position position)
			throw();
	
		/** Get a constant pointer to a child SecondaryStructure at a given position.
				The pointer is 0 if this instance does not have a SecondaryStructure at this position.
				@param   position the position of the child SecondaryStructure
				@return  SecondaryStructure* -
								 constant pointer to the child SecondaryStructure at <b>  position </b>
		*/
		const SecondaryStructure* getSecondaryStructure(Position position) const
			throw();
	
		/** Get a pointer to a child Residue at a given position.
				The pointer is 0 if this instance does not have a Residue at this position.
				@param   position the position of the child Residue
				@return  Residue* - mutable pointer to the child Residue at <b>  position </b>
		*/
		Residue* getResidue(Position position)
			throw();
	
		/** Get a pointer to a child Residue at a given position.
				The pointer is 0 if this instance does not have a Residue at this position.
				@param   position the position of the child Residue
				@return  Residue* - constant pointer to the child Residue at <b>  position </b>
		*/
		const Residue* getResidue(Position position) const
			throw();
	
		/** Get a pointer to the N-terminal Residue.
				The pointer is 0 if this instance does not have a Residue with
				the property "PROPERTY__AMINO_ACID".
				@return  Residue* - mutable pointer to the N-terminal Residue
		*/
		Residue* getNTerminal()
			throw();
	
		/** Get a constant pointer to the N-terminal Residue.
				The pointer is 0 if this instance does not have a Residue with
				the property "PROPERTY__AMINO_ACID".
				@return  Residue* - constant pointer to the N-terminal Residue
		*/
		const Residue* getNTerminal() const
			throw();

		/** Get a pointer to the C-terminal Residue.
				The pointer is 0 if this instance does not have a Residue with
				the property "PROPERTY__AMINO_ACID".
				@return  Residue* - mutable pointer to the C-terminal Residue
		*/
		Residue* getCTerminal()
			throw();
	
		/** Get a constant pointer to the C-terminal Residue.
				The pointer is 0 if this instance does not have a Residue with
				the property "PROPERTY__AMINO_ACID".
				@return  Residue* - constant pointer to the C-terminal Residue
		*/
		const Residue* getCTerminal() const
			throw();

		/** Get a pointer to a child PDBAtom at a given position.
				The pointer is 0 if this instance does not have a PDBAtom at this position.
				@param   position the position of the child PDBAtom
				@return  PDBAtom* - mutable pointer to the child PDBAtom at <b>  position </b>
		*/
		PDBAtom* getPDBAtom(Position position)
			throw();
	
		/** Get a pointer to a child PDBAtom at a given position.
				The pointer is 0 if this instance does not have a PDBAtom at this position.
				@param   position the position of the child PDBAtom
				@return  PDBAtom* - constant pointer to the child PDBAtom at <b>  position </b>
		*/
		const PDBAtom* getPDBAtom(Position position) const
			throw();
	
		/** Count the SecondaryStructures
				@return  Size the number of secondary structures
		*/
		Size countSecondaryStructures() const
			throw();

		/** Count the Residues
				@return  Size the number of residues
		*/
		Size countResidues() const
			throw();

		/** Count the PDBAtoms
				@return  Size the number of PDBAtoms
		*/
		Size countPDBAtoms() const
			throw();

		/** Prepend a SecondaryStructure at position 0.
				@param secondary_structure the SecondaryStructure to prepend
		*/
		void prepend(SecondaryStructure& secondary_structure)
			throw();

		/** Append a SecondaryStructure after the last position.
				@param secondary_structure the SecondaryStructure to append
		*/
		void append(SecondaryStructure& secondary_structure)
			throw();

		/** Insert a SecondaryStructure after the last position.
				@param secondary_structure the SecondaryStructure to insert
		*/
		void insert(SecondaryStructure& secondary_structure)
			throw();

		/** Insert a SecondaryStructure before a given Composite object.
				@param secondary_structure the SecondaryStructure to insert
				@param before the Composite object to insert before
		*/
		void insertBefore(SecondaryStructure& secondary_structure, Composite& before)
			throw();

		/** Insert a SecondaryStructure after a given Composite object.
				@param secondary_structure the SecondaryStructure to insert
				@param after the Composite object to insert after
		*/
		void insertAfter(SecondaryStructure& secondary_structure, Composite& after)
			throw();

		/** Remove a SecondaryStructure
				@param secondary_structure the SecondaryStructure to remove
		*/
		bool remove(SecondaryStructure& secondary_structure)
			throw();

		/** Prepend a Residue at position 0.
				@param residue the Residue to prepend
		*/
		void prepend(Residue& residue)
			throw();

		/** Append a Residue after the last position.
				@param residue the Residue to append
		*/
		void append(Residue& residue)
			throw();

		/** Insert a Residue after the last position.
				@param residue the Residue to insert
		*/
		void insert(Residue& residue)
			throw();

		/** Insert a Residue before a given Composite object.
				@param residue the Residue to insert
				@param before the Composite object to insert before
		*/
		void insertBefore(Residue& residue, Composite& before)
			throw();

		/** Insert a Residue after a given Composite object.
				@param residue the Residue to insert
				@param after the Composite object to insert after
		*/
		void insertAfter(Residue& residue, Composite& after)
			throw();

		/** Remove a Residue.
				@param residue the Residue to remove
		*/
		bool remove(Residue& residue)
			throw();

		/**	Cut all children of <tt>chain</tt> and prepend them before the children of this chain.
				@param chain the chain to access
		*/
		void spliceBefore(Chain& chain)
			throw();

		/**	Cut all children of <tt>chain</tt> and append them after the children of this chain.
				@param chain the chain to access
		*/
		void spliceAfter(Chain &chain)
			throw();

		/**	Move the children of <tt>chain</tt> into this chain.
				The children of <tt>chain</tt> are inserted using  \link spliceBefore spliceBefore \endlink .
		*/
		void splice(Chain &chain)
			throw();
		
		//@}
		/**	@name	Debugging and Diagnostics 
		*/
		//@{
		
		/** Internal state dump.
				Dump the current internal state to the output 
				ostream <b>  s </b> with dumping depth <b>  depth </b>.
				@param   s - output stream where to output the internal state.
				@param   depth - the dumping depth
		*/
		virtual void dump(std::ostream& s = std::cout, Size depth = 0) const
			throw();

		//@}

		// --- EXTERNAL ITERATORS

		BALL_KERNEL_DEFINE_ITERATOR_CREATORS(SecondaryStructure)
		BALL_KERNEL_DEFINE_ITERATOR_CREATORS(Residue)
		BALL_KERNEL_DEFINE_ITERATOR_CREATORS(PDBAtom)

	protected:

	private:

		AtomContainer* getAtomContainer(Position position)
			throw();
	
		const AtomContainer* getAtomContainer(Position position) const
			throw();
	
		Atom* getAtom(Position position)
			throw();
	
		const Atom* getAtom(Position position) const
			throw();

		void prepend(Atom& atom)
			throw();

		void append(Atom& atom)
			throw();

		void insert(Atom& atom)
			throw();

		void insertBefore(Atom& atom, Composite& before)
			throw();

		void insertAfter(Atom& atom, Composite& after)
			throw();

		bool remove(Atom& atom)
			throw();

		void prepend(AtomContainer& atom_container)
			throw();

		void append(AtomContainer& atom_container)
			throw();

		void insert(AtomContainer& atom_container)
			throw();

		void insertBefore(AtomContainer& atom_container, Composite& before)
			throw();

		void insertAfter(AtomContainer& atom_container, Composite& after)
			throw();

		void spliceBefore(AtomContainer& atom_container)
			throw();

		void spliceAfter(AtomContainer& atom_container)
			throw();

		void splice(AtomContainer& atom_container)
			throw();

		bool remove(AtomContainer& atom_container)
			throw();

		BALL_KERNEL_DEFINE_ITERATOR_CREATORS(AtomContainer)
	};
} // namespace BALL

#endif // BALL_KERNEL_CHAIN_H
