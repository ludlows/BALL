// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// $Id: claverieParameter.h,v 1.14.6.1 2005/07/29 12:38:01 amoll Exp $
//

#ifndef BALL_SOLVATION_CLAVERIEPARAMETER_H
#define BALL_SOLVATION_CLAVERIEPARAMETER_H

#ifndef BALL_FORMAT_PARAMETERS_H
# include <BALL/FORMAT/parameters.h>
#endif

#ifndef BALL_FORMAT_PARAMETERSECTION_H
# include <BALL/FORMAT/parameterSection.h>
#endif

#ifndef BALL_MOLMEC_PARAMETER_FORCEFIELDPARAMETERS_H
# include <BALL/MOLMEC/PARAMETER/forceFieldParameters.h>
#endif

namespace BALL
{
	/** ClaverieParameter class. 
			This class provides the parameters needed for the computation of
			van-der-Waals interaction energies according to the scheme by Huron and
			Claverie. [missing: units, reference].
			<b>Note:</b> Dispersion and repulsion terms are <b>not</b> distinguished.  \par
		\ingroup Solvation	
	 */
	class BALL_EXPORT ClaverieParameter
		:	public ParameterSection
	{

		public:

		BALL_CREATE(ClaverieParameter)

		/** @name Constructors and destructors. 
		*/
		//@{

		/** Default constructor 
		*/
		ClaverieParameter() throw();

		/** Detailed constructor 
		*/
		ClaverieParameter(const Parameters& parameters) 
			throw();

		/** Copy constructor 
		*/
		ClaverieParameter(const ClaverieParameter& param) 
			throw();

		/** Destructor 
		*/
		virtual ~ClaverieParameter() 
			throw();

		//@}
		/** @name Assignment 
		*/
		//@{

		/** Assignment operator 
		*/
		const ClaverieParameter& operator = (const ClaverieParameter& param)
			throw();

		/** Clear method 
		*/
		virtual void clear() 
			throw();

		//@}
		/** @name Accessors 
		*/
		//@{

		/** Indicate, whether this parameter set has parameters for the
				specified atom types.
		*/
		bool hasParameters(Atom::Type solvent_type, Atom::Type solute_type) const
			throw();

		/** Get the parameters for the specified atom types.
		*/
		std::pair<float, float> getParameters(Atom::Type solvent_type,
				Atom::Type solute_type) const 
			throw();

		/** ?????
		*/
		std::pair<float, float> getParameters(Atom::Type type) const 
			throw();

		/** Return the whole parameter set.
		*/
		const ::std::vector< ::std::pair<float, float> >& getParameters() const
			throw();

		/** Return the indices used for mapping types to numbers
		*/
		const HashMap<Atom::Type, Index>& getIndices() const
			throw();

		//@}
		/** @name Predicates 
		*/
		//@{

		/** Equality operator.
		*/
		bool operator == (const ClaverieParameter& param) const 
			throw();

		//@}

		/** Extract the parameter file section (@see ParameterSection).
		*/
		virtual bool extractSection(ForceFieldParameters& parameters,
				const String& section_name) 
			throw();

		protected:

		/*_ This vector contains the paramaters that were read from the
			  parameter file 
		*/
		std::vector< std::pair<float, float> > parameters_;

		/*_ Here the atom types are mapped to the indices of the vector 
		*/
		HashMap<Atom::Type, Index> indices_;

	};
   
} // namespace BALL

#endif // BALL_SOLVATION_CLAVERIEPARAMETER_H
