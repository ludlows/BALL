// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// $Id: forceFieldParameters.h,v 1.19.6.1 2005/07/29 12:37:49 amoll Exp $
//

// Molecular Mechanics: general force field parameter class

#ifndef BALL_MOLMEC_FORCEFIELDPARAMETERS_H
#define BALL_MOLMEC_FORCEFIELDPARAMETERS_H

#ifndef BALL_COMMON_H
#	include <BALL/common.h>
#endif

#ifndef BALL_FORMAT_PARAMETERS_H
# include <BALL/FORMAT/parameters.h>
#endif

#ifndef BALL_MOLMEC_PARAMETER_ATOMTYPES_H
# include <BALL/MOLMEC/PARAMETER/atomTypes.h>
#endif

namespace BALL 
{
	class AtomTypes;
	
	/**	Force field parameter class.
			
			\ingroup  MolmecParameters
	*/
	class BALL_EXPORT ForceFieldParameters
		:	public Parameters
	{
		public:

		BALL_CREATE(ForceFieldParameters)

		friend class ForceField;

		/**@name	Constructors and destructor	
		*/
		//@{

		/**	Default constructor.
		*/
		ForceFieldParameters();

		/**	Constructor.
		*/
		ForceFieldParameters(const String& filename);

		/**	Copy constructor
		*/
		ForceFieldParameters(const ForceFieldParameters& force_field_parameter);

		/**	Destructor.
		*/
		virtual ~ForceFieldParameters()
			throw();

		//@}
		/** @name Assignment 
		*/
		//@{

		/** Clear method 
		*/
		virtual void clear()
			throw();

		/** Assignment operator 
		*/
		const ForceFieldParameters& operator = (const ForceFieldParameters& param);
		
		//@}
		/**@name	Accessors 	
		*/
		//@{

		/**	Return a reference to the atom type parameter section
		*/
		AtomTypes&	getAtomTypes();

		/**	Read the contents of the INI file and interpret them.
		*/
		virtual bool init();

		//@}
		/**	@name	Predicates
		*/
		//@{
			
		/**	Valididty predicate.
				Return <b>true</b> if the force field parameters were correctly
				initialized, the internal INI file is valid and the internal atom types		
				object is valid.
				@return bool - <tt>valid_ && parameter_file_.isValid() && atom_types_.isValid()</tt>
		*/
		virtual bool isValid() const;

		/** Equality operator 
		*/
		bool operator == (const ForceFieldParameters& param) const;

		//@}

		protected:

		/*_@name	Protected Members 
		*/
		//_@{ 

		/*_	the atom types section
		*/
		AtomTypes	atom_types_;

		//_@} 
	};
} // namespace BALL

#endif // BALL_MOLMEC_FORCEFIELDPARAMETERS_H
