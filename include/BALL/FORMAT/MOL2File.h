// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// $Id: MOL2File.h,v 1.22.6.2 2005/08/12 12:49:43 amoll Exp $
//

#ifndef BALL_FORMAT_MOL2FILE_H
#define BALL_FORMAT_MOL2FILE_H

#ifndef BALL_FORMAT_GENERICMOLFILE_H
#	include <BALL/FORMAT/genericMolFile.h>
#endif

#ifndef BALL_MATHS_VECTOR3_H
#	include <BALL/MATHS/vector3.h>
#endif

namespace BALL 
{
	class Atom;
	class System;

	/**	SYBYL MOL2 file class.
			This class is used to read and write SYBYL MOL2 files (Tripos).	 \par
			
    	\ingroup  StructureFormats
	*/
	class BALL_EXPORT MOL2File
		: public GenericMolFile
	{
		public:

		/**	@name	Constants	
		*/
		//@{
			
		/**	Tripos Record Type Identifier (RTI)
		*/
		static const String TRIPOS;
		//@}
		
		/**	@name	Constructors and Destructors
		*/
		//@{

		/**	Default constructor
		*/
		MOL2File()
			throw();

		/** Detailed constructor
		*/
		MOL2File(const String& filename, File::OpenMode open_mode = std::ios::in)
			throw(Exception::FileNotFound);

		/**	Copy constructor
		*/
		MOL2File(const MOL2File& file)
			throw(Exception::FileNotFound);

		/// Destructor
		virtual ~MOL2File()
			throw();
		
		//@}

		/**	@name Reading and Writing of Kernel Datastructures
		*/
		//@{
		
		/**	Write a system to the MOL2 file
		*/
		virtual bool write(const System& system)
			throw(File::CannotWrite);
		
		/**	Read a system from the MOL2 file
		*/
		virtual bool read(System&	system)
			throw(Exception::ParseError);

		/**	
		*/
		virtual Molecule* read()
			throw(Exception::ParseError);

		/**	
		*/
		virtual bool write(const Molecule& molecule)
			throw(File::CannotWrite);

		///
		const MOL2File& operator = (const MOL2File& file) throw();

		//@}

		protected:
		
		void readAtomSection_();

		void readBondSection_();

		void  readMoleculeSection_();

		void readSetSection_();

		void readSubstructureSection_();

		String getSybylType_(const Atom& atom) const;

		bool nextLine_();
		
		void clear_();
		
		bool buildAll_(System& system);

		BALL_EXPORT struct AtomStruct
		{
			String		name;
			Vector3		position;
			String		type;
			Position	substructure;
			String		substructure_name;
			float			charge;
		};

		BALL_EXPORT struct BondStruct
		{
			Position	atom1;
			Position	atom2;
			String		type;
		};

		BALL_EXPORT struct MoleculeStruct
		{
			String			name;
			Size				number_of_atoms;
			Size				number_of_bonds;
			Size				number_of_substructures;
			Size				number_of_features;
			Size				number_of_sets;
			String			type;
			String			charge_type;
			String			comment;
		};

		BALL_EXPORT struct SetStruct
		{
			String			name;
			String			type;
			String			subtype;
			String			comment;
			Size				number_of_members;
			vector<int>	members;
		};
		
		BALL_EXPORT struct SubstructureStruct
		{
			String			name;
			Size				root_atom;
			String			substructure_type;
			Size				dictionary_type;
			String			chain;
			String			sub_type;
			Size				inter_bonds;
			String			comment;
		};

			
		vector<AtomStruct>					atoms_;
		vector<BondStruct>					bonds_;
		vector<SetStruct>						sets_;
		vector<SubstructureStruct>	substructures_;
		MoleculeStruct							molecule_;

		Size	number_of_lines_;
		static const Size MAX_LENGTH_;
		char buffer_[4096];
		String line_;
	};
} // namespace BALL

#endif // BALL_FORMAT_MOL2FILE_H
