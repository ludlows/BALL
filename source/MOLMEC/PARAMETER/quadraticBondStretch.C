// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// $Id: quadraticBondStretch.C,v 1.10 2003/08/26 09:17:54 oliver Exp $
//

#include <BALL/MOLMEC/PARAMETER/quadraticBondStretch.h>
#include <BALL/MOLMEC/PARAMETER/forceFieldParameters.h>

using namespace std;

namespace BALL 
{

	QuadraticBondStretch::QuadraticBondStretch()
		:	ParameterSection(),
			k_(0),
			r0_(0),
			is_defined_(0)
	{
	}

	QuadraticBondStretch::~QuadraticBondStretch() throw()
	{
		clear();
	}

	void QuadraticBondStretch::clear() throw()
	{
		// clear allocated parameter fields
		delete [] k_;
		delete [] r0_;
		delete [] is_defined_;

		k_ = 0;
		r0_ = 0;
		is_defined_ = 0;

		ParameterSection::clear();
	}

	bool QuadraticBondStretch::extractSection
		(Parameters& parameters, const String& section_name)
	{
		return ParameterSection::extractSection(parameters, section_name);
	}

	bool QuadraticBondStretch::extractSection
		(ForceFieldParameters& parameters, const String& section_name)
	{
		// check whether the parameters are valid
		if (!parameters.isValid())
		{
			return false;
		}
			

		// extract the basis information
		ParameterSection::extractSection(parameters, section_name);
		
		// check whether all variables we need are defined, terminate otherwise
		if (!hasVariable("r0") || !hasVariable("k"))
		{
			return false;
		}

		// build a two dimensional array of the atom types
		// loop variable
		Size	i;

		AtomTypes& atom_types = parameters.getAtomTypes();
		number_of_atom_types_ = atom_types.getNumberOfTypes();
		
		// allocate two onedimensional fields for the two parameters
		// and a two dimensional field of boolean variables.
		// we might have to delete old stuff lying around
		if (k_)
			delete [] k_;
		if (r0_)
			delete [] r0_;
		if (is_defined_)
			delete [] is_defined_;

		k_  = new float[number_of_atom_types_ * number_of_atom_types_];
		r0_ = new float[number_of_atom_types_ * number_of_atom_types_];
		is_defined_ = new bool[number_of_atom_types_ * number_of_atom_types_];

		for (i = 0; i < number_of_atom_types_ * number_of_atom_types_; i++) 
		{
			is_defined_[i] = false;
		}

		StringHashMap<Index>::Iterator it;

		// determine the factor to convert the parameters to the standard units used
		// as a default, energies are assumend to be in kJ/mol and distances in Angstrom
		double factor_k = 1.0;
		double factor_r0 = 1.0;
		if (options.has("unit_k"))
		{
			if (options["unit_k"] == "kcal/mol")
			{
				factor_k = Constants::JOULE_PER_CAL;
			}

			if (options["unit_k"] == "cal/mol")
			{
				factor_k = Constants::JOULE_PER_CAL * 0.001;
			}
			
			if (options["unit_k"] == "J/mol")
			{
				factor_k = 0.001;
			}			
		}	
		
		if (options.has("unit_r0"))
		{
			if (options["unit_r0"] == "pm")
			{
				factor_k = 0.1;
			}
		}	
		

		Atom::Type		type_I;
		Atom::Type		type_J;
		String				type_name_I;
		String				type_name_J;
		String				key;
		Index					index = 0;

		for (it = section_entries_.begin(); !(it == section_entries_.end()); ++it)
		{
			key = (*it).first;
			if ((key.size() > 0) && (key.find_first_of(" ", 0) > 0)) 
			{
				type_name_I = key.before(" ", 0);
				type_name_J = key.after(" ", 0);
				if ((atom_types.hasType(type_name_I)) && (atom_types.hasType(type_name_J))) 
				{
					type_I = atom_types.getType(type_name_I);
					type_J = atom_types.getType(type_name_J);
					index = (Index)(type_I * number_of_atom_types_ + type_J);
					is_defined_[index] = true;
					k_ [index] = getValue(key, "k").toFloat() * factor_k;
					r0_ [index] = getValue(key, "r0").toFloat() * factor_r0;
					index = (Index)(type_I + number_of_atom_types_ * type_J);
					is_defined_[index] = true;
					k_ [index] = getValue(key, "k").toFloat() * factor_k;
					r0_ [index] = getValue(key, "r0").toFloat() * factor_r0;
				}
			}
		}

		return true;
	}


	bool QuadraticBondStretch::hasParameters(Atom::Type I, Atom::Type J) const 
	{
		if ((I < 0) || ((Size)I >= number_of_atom_types_))
		{
			return false;
		}

		if ((J < 0) || ((Size)J >= number_of_atom_types_))
		{
			return false;
		}

		return is_defined_[I * number_of_atom_types_ + J];
	}


	QuadraticBondStretch::Values QuadraticBondStretch::getParameters
		(Atom::Type I, Atom::Type J) const 
	{
		QuadraticBondStretch::Values parameters;
		assignParameters(parameters, I, J);
		return parameters;
	}


	bool QuadraticBondStretch::assignParameters
		(QuadraticBondStretch::Values& parameters,
		 Atom::Type I, Atom::Type J) const 
	{
		if (hasParameters(I, J)) 
		{
			parameters.k = k_[I * number_of_atom_types_ + J];
			parameters.r0 = r0_[I * number_of_atom_types_ + J];

			return true;
		}

		return false;
	}
	 
} // namespace BALL
