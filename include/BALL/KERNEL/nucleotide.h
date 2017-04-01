// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// $Id: nucleotide.h,v 1.36.4.2 2005/08/11 14:51:58 amoll Exp $
//

#ifndef BALL_KERNEL_NUCLEOTIDE_H
#define BALL_KERNEL_NUCLEOTIDE_H

#ifndef BALL_KERNEL_FRAGMENT_H
#	include <BALL/KERNEL/fragment.h>
#endif

#ifndef BALL_KERNEL_RESIDUE_H
#	include <BALL/KERNEL/residue.h>
#endif

#ifndef BALL_KERNEL_NUCLEOTIDEITERATOR_H
#	include <BALL/KERNEL/nucleotideIterator.h>
#endif

 
#define BALL_NUCLEOTIDE_DEFAULT_ID               ""
#define BALL_NUCLEOTIDE_DEFAULT_INSERTION_CODE   ' '

namespace BALL 
{
	class NucleicAcid;

	/**	Nucleotide class.
			This class is used to represent nucleotides within
			an  \link NucleicAcid NucleicAcid \endlink  object.
			 \par
			
    	\ingroup KernelContainers 
	*/
	class BALL_EXPORT Nucleotide
		: public Fragment
	{
		public:

		BALL_CREATE_DEEP(Nucleotide)
	
		/**	@name	Enums
		*/
		//@{

		///
		enum Property
		{
			///
			PROPERTY__5_PRIME = Residue::NUMBER_OF_PROPERTIES + 1,
			///
			PROPERTY__3_PRIME,
			/// 
			PROPERTY__NUCLEOTIDE,

			///
			NUMBER_OF_PROPERTIES
		};

		//@}
		/**	@name	Constructors and Destructors 
		*/
		//@{	
	
		/// Default constructor
		Nucleotide()
			throw();
	
		/// Copy constructor
		Nucleotide(const Nucleotide& nucleotide, bool deep = true)
			throw();
	
		/// Detailled constructor
		Nucleotide
			(const String& name,
			 const String& id = BALL_NUCLEOTIDE_DEFAULT_ID,
			 char insertion_code = BALL_NUCLEOTIDE_DEFAULT_INSERTION_CODE)
			throw();

		/// Destructor
		virtual ~Nucleotide()
			throw();

		/// Clears the nucleotides contents.
		virtual void clear()
			throw();
	
		/// Clears the nucleotides contents and removes it from all composite structures.
		virtual void destroy()
			throw();
	
		//@}
		/** @name Persistence 
		*/
		//@{

		/**	Writes a Nucleotide object to a persistent stream.
				@param pm the persistence manager
		*/
		void persistentWrite(PersistenceManager& pm, const char* name = 0) const
			throw(Exception::GeneralException);

		/**	Reads a Nucleotide object from a persistent stream.
				@param pm the persistence manager
		*/
		void persistentRead(PersistenceManager& pm)
			throw(Exception::GeneralException);

		//@}
		/**	@name	Assignment 
		*/
		//@{

		/** Assignment with cloning facility.
				The assignment is either deep or shallow (default).
				@param  nucleotide the nucleotide to be copied (cloned)
				@param  deep make a deep (=<tt>true</tt>) or shallow (=<tt>false</tt>) copy
		*/
		void set(const Nucleotide& nucleotide, bool deep = true)
			throw();

		/** Assignment operator.
				The assignment is either deep or shallow (default).
				@param   nucleotide the nucleotide to be copied (cloned)
				@return  nucleotide& - this instance nucleotide
				@see     nucleotide::set
		*/
		Nucleotide& operator = (const Nucleotide& nucleotide)
			throw();

		/** Copy this instance to <b>  nucleotide </b>.
				The assignment is either deep or shallow (default).
				@param  nucleotide the nucleotide to be assigned to
				@see    nucleotide::set
		*/
		void get(Nucleotide& nucleotide, bool deep = true) const
			throw();

		/** Swapping of instaces of nucleotide.
				@param  nucleotide the instance of nucleotide to swap with
		*/
		void swap(Nucleotide& nucleotide)
			throw();
	
		//@}

		/**	Equality operator.
				Two nucleotides are equal if they have the same handle.
				@see Object::operator ==.
		*/
		bool operator == (const Nucleotide& nucleotide) const
			throw();

		/**	Inequality operator
				@see operator ==
		*/
		bool operator != (const Nucleotide& nucleotide) const
			throw();

		/**	@name	Accessors */
		//@{

		/** Get a pointer to the parent NucleicAcid.
				The pointer is 0 if this instance nucleotide does not have a parent NucleicAcid.
				@return  NucleicAcid* - mutable pointer to the parent NucleicAcid
		*/
		NucleicAcid* getNucleicAcid()
			throw();
		
		/** Get a pointer to the parent NucleicAcid.
				The pointer is 0 if this instance nucleotide does not have a parent NucleicAcid.
				@return  NucleicAcid* - constant pointer to the parent NucleicAcid
		*/
		const NucleicAcid* getNucleicAcid() const
			throw();

		/**	Set the ID of the nucleotide.
				@param id the new ID
		*/
		void setID(const String& id)
			throw();

		/**	Retrieve the ID of the nucleotide.
				@return String the ID (constant)
		*/
		const String& getID() const
			throw();

		/**	Set the insertion code of the nucleotide.
				@param insertion_code the new insertion code
		*/
		void setInsertionCode(char insertion_code)
			throw();

		/**	Retrieve the insertion code of the nucleotide.
				@return String the insertion code (constant)
		*/
		char getInsertionCode() const
			throw();

		/** Prepend an atom at position 0.
				@param atom the atom to prepend
		*/
		void prepend(Atom& atom)
			throw();

		/** Append an atom after the last position.
				@param atom the atom to append
		*/
		void append(Atom& atom)
			throw();

		/** Insert an atom after the last position.
				@param atom the atom to insert
		*/
		void insert(Atom& atom)
			throw();

		/** Insert an atom before a given Composite object.
				@param atom the atom to insert
				@param before the Composite object to insert before
		*/
		void insertBefore(Atom& atom, Composite& before)
			throw();

		/** Insert an atom after a given Composite object.
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

		/**	Cut all children of <tt>nucleotide</tt> and prepend them before the children of this instance.
				@param nucleotide the nucleotide to access
		*/
		void spliceBefore(Nucleotide& nucleotide)
			throw();

		/**	Cut all children of <tt>nucleotide</tt> and append them after the children of this instance.
				@param nucleotide the nucleotide to access
		*/
		void spliceAfter(Nucleotide& nucleotide)
			throw();

		/**	Move the children of <tt>nucleotide</tt> into this instance.
				The children are inserted using  \link spliceBefore spliceBefore \endlink .
		*/
		void splice(Nucleotide& nucleotide)
			throw();

		//@}
		/**	@name	Predicates 
		*/
		//@{

		/**	Test if this instance nucleotide is terminal.
				Returns true, if this instance is the first or 
				last nucleotide in its parent NucleicAcid.
				@return bool
		*/
		bool isTerminal() const
			throw();
	
		/**	Test if this instance nucleotide is 3-prime.
				Returns true, if this instance is the first nucleotide in its parent NucleicAcid.
				@return bool
		*/
		bool is3Prime() const
			throw();

		/**	Test if this instance nucleotide is 5-prime.
				Returns true, if this instance is the last nucleotide in its parent NucleicAcid.
				@return bool
		*/
		bool is5Prime() const
			throw();

		//@}
		/**	@name	Debugging and Diagnostics 
		*/
		//@{

		/** Internal state and consistency self-validation.
				@return	 bool -	<tt>true</tt> if the internal state of this instance nucleotide is correct
												(self-validated) and consistent, <tt>false</tt> otherwise
		*/
		virtual bool isValid() const
			throw();

		/** Internal state dump.
				Dump the current internal state of this instance 
				to the output ostream <b>  s </b> with dumping depth <b>  depth </b>.
				@param	s output stream where to output the internal state of this instance nucleotide
				@param  depth the dumping depth
		*/
		virtual void dump(std::ostream& s = std::cout, Size depth = 0) const
			throw();

		//@}
 
		private:

		AtomContainer* getAtomContainer(Position position)
			throw();
	
		const AtomContainer* getAtomContainer(Position position) const
			throw();
	
		Size countAtomContainers() const
			throw();

		void prepend(AtomContainer& atom_container)
			throw();

		void append(AtomContainer& atom_container)
			throw();

		void insert(AtomContainer& atom_container)
			throw();

		void insertBefore(AtomContainer& atom_container, Composite& composite)
			throw();

		void insertAfter(AtomContainer& atom_container, Composite& composite)
			throw();

		void spliceBefore(AtomContainer& atom_container)
			throw();

		void spliceAfter(AtomContainer& base_ragment)
			throw();

		void splice(AtomContainer& AtomContainer)
			throw();

		bool remove(AtomContainer& AtomContainer)
			throw();

		bool isSuperAtomContainerOf(const AtomContainer& atom_container) const
			throw();

		BALL_KERNEL_DEFINE_ITERATOR_CREATORS(AtomContainer)

		// --- ATTRIBUTES

		String 	id_;

		char 		insertion_code_;
	};


  template <class NucleotideContainerType>
  const Nucleotide* get5Prime(const NucleotideContainerType& nucleotide_container)
  {
		NucleotideConstIterator res_it;
    for ( res_it = nucleotide_container.beginNucleotide(); !res_it.isEnd(); ++res_it)
		{
			return &(*res_it);
		}

    return 0;
  }

  template <class NucleotideContainerType>
  const Nucleotide* get3Prime(const NucleotideContainerType& nucleotide_container)
  {
    for (NucleotideConstIterator res_it = nucleotide_container.rbeginNucleotide(); !res_it.isREnd(); ++res_it)
		{
			return &(*res_it);
		}

    return 0;
  }
 
} // namespace BALL

#endif // BALL_KERNEL_NUCLEOTIDE_H
