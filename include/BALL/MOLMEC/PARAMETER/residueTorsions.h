// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// $Id: residueTorsions.h,v 1.19.4.2 2005/08/11 14:15:55 amoll Exp $
//

// Molecular Mechanics Parameter: class describing the ResidueTorsions section of a parameter file
 
#ifndef BALL_MOLMEC_PARAMETER_RESIDUETORSIONS_H
#define BALL_MOLMEC_PARAMETER_RESIDUETORSIONS_H

#ifndef BALL_FORMAT_PARAMETERSECTION_H
#	include <BALL/FORMAT/parameterSection.h>
#endif

#ifndef BALL_MOLMEC_PARAMETER_ATOMTYPES_H
#	include <BALL/MOLMEC/PARAMETER/atomTypes.h>
#endif

namespace BALL 
{
	/**	Parameter class containing all proper torsions occuring in a residue.
			Several force fields (e.g. CHARMM) do not necessarily consider or
			parametrize all occurring torsions but explicitly list the torsions
			for each residue. This parameter section is used to represent this list.
			The class AMBER and CHARMM torsions components check for the presence of
			the parameter section [ResidueTorsions] and then decide whether they have
			to generate the torsions by itself (creating all torsions and
			complaining about missing parameters) or whether they have to read them
			from this section. \par
			
    	\ingroup  MolmecParameters
	*/
	class BALL_EXPORT ResidueTorsions 
		:	public ParameterSection
	{
		public:

		/**	@name	Type definitions
		*/
		//@{

		/**	Strcuture containing the names of the residue and the atoms for a torsion.	
		*/
		struct BALL_EXPORT Data
		{
			String	residue_name;
			String	atom_name_A;
			String	atom_name_B;
			String	atom_name_C;
			String	atom_name_D;
			
			Data(const String& name, const String& A, const String& B, const String& C, const String& D)
				:	residue_name(name),
					atom_name_A(A),
					atom_name_B(B),
					atom_name_C(C),
					atom_name_D(D)
			{
			}

			Data()
				:	residue_name(""),
					atom_name_A(""),
					atom_name_B(""),
					atom_name_C(""),
					atom_name_D("")
			{
			}

			bool operator == (const Data& data) const;
			bool operator != (const Data& data) const;
		};

		//@}
		/**	@name	Constructors and Destructors
		*/
		//@{

		/**	Default constructor.
		*/
		ResidueTorsions();

		/**	Destructor.
		*/
		virtual ~ResidueTorsions() throw();

		/**	Clear method.  
		*/
		virtual void clear() throw();

		//@}
		/**	@name	Accessors
		*/
		//@{

		/**	Reads a parameter section from an INI file.
				This method reads the section given in section_name from ini_file,
				interprets (if given) a format line, reads the data from this section according to 
				the format, and builds some datastructures for fast and easy acces this data.
		*/
		virtual bool extractSection(ForceFieldParameters& parameters, const String& section_name);

		///
		virtual bool extractSection(Parameters& parameters, const String& section_name);

		/**	Return the number of torsions for this residue.
		*/
		Size getNumberOfResidueTorsions(const String& residue_name) const;

		/**	Assign the <i>  i </i>th torsion for a residue.
				@param	name the residue name (including modifiers like -S or -N);
				@param  i the index. 0 $<$ <tt>i</tt> $<$  \link getNumberOfResidueTorsions getNumberOfResidueTorsions \endlink 
				@param	ResidueTorsion the torsion to be assigned to
				@return bool - <b>true</b> if the torsion was found, <b>false</b> otherwise
		*/
		bool assignTorsion(const String& name, Position i, Data& torsion) const;

		/**	Return true if the torsion has to be considered for the residue.
		*/
		bool hasTorsion
			(const String& residue, const String& atom_A, const String& atom_B,
			 const String& atom_C, const String& atom_D) const;

		//@}

		protected:

		/*_	Contains arrays of ResidueTorsions.
				All torsions for a given residue name are collected in 
				a vector and accessed via the residue name through a StringHashMap.
		*/
		StringHashMap<vector<Data> >	torsions_;

		/*_ Hash set of all torsion identifiers.
				This hash set contains all entries in the form of strings.
				It is used by \Ref{hasTorsion}.
		*/
		HashSet<String>								all_torsions_;
	};
} // namespace BALL

#endif // BALL_MOLMEC_PARAMETER_RESIDUETORSIONS_H
