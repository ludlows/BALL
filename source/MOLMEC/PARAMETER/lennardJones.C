// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// $Id: lennardJones.C,v 1.19 2005/02/18 13:06:53 amoll Exp $
//

#include <BALL/MOLMEC/PARAMETER/lennardJones.h>
#include <BALL/MOLMEC/PARAMETER/forceFieldParameters.h>

using namespace std;

namespace BALL 
{

	LennardJones::LennardJones() throw()
		:	ParameterSection(),
			A_(0),
			B_(0),
			N_(0),
			Aij_(0),
			Bij_(0),
			is_defined_(0),
			format_(EPSILON_R_FORMAT),
			names_()
	{
	}


	LennardJones::LennardJones(const LennardJones& lj) throw()
		:	ParameterSection(lj),
			A_(lj.A_),
			B_(lj.B_),
			N_(lj.N_),
			Aij_(lj.Aij_),
			Bij_(lj.Bij_),
			is_defined_(lj.is_defined_),
			format_(lj.format_),
			names_(lj.names_)
	{
	}


	LennardJones::~LennardJones() throw()
	{
		clear();
	}


	void LennardJones::clear() throw()
	{
		// clear allocated parameter fields
		A_.clear();
		B_.clear();
		N_.clear();
		Aij_.clear();
		Bij_.clear();
		is_defined_.clear();
		format_ = EPSILON_R_FORMAT;
		names_.clear();

		ParameterSection::clear();
	}


	const LennardJones& LennardJones::operator = (const LennardJones& lj)
		throw()
	{
		ParameterSection::operator = (lj);
		A_ = lj.A_;
		B_ = lj.B_;
		N_ = lj.N_;
		Aij_ = lj.Aij_;
		Bij_ = lj.Bij_;
		is_defined_ = lj.is_defined_;
		format_ = lj.format_;
		names_ = lj.names_;

		return *this;
	}


	bool LennardJones::extractSection
		(Parameters& parameters, const String& section_name) throw()
	{
		return ParameterSection::extractSection(parameters, section_name);
	}


	bool LennardJones::extractSection
		(ForceFieldParameters& parameters, const String& section_name) throw()
	{
		// check whether the parameters are valid
		if (!parameters.isValid())
		{
			return false;
		}
		
		// extract the basis information
		ParameterSection::extractSection(parameters, section_name);
		bool use_geometric_mean = false;


		// check whether all variables we need are defined, terminate otherwise
		if ((!hasVariable("A") || !hasVariable("B"))
				&& (!hasVariable("epsilon") || !hasVariable("R")))
		{
			Log.error() << "LennardJones::extractSection: Lennard Jones parameter section requires two variable columns:"		
				<< "A/B or epsilon/R" << endl;

			return false;

		} 
		else 
		{
			// format_ == A_B_FORMAT:				parameters are in A/B format
			// format_ == EPSILON_R_FORMAT:	parameters are in epsilon/R format
			// format_ == SLATER_KIRKWOOD_FORMAT:	parameters are in epsilon/R format
			if (hasVariable("epsilon") && hasVariable("R"))
			{
				format_ = EPSILON_R_FORMAT;
				if (options.has("radius_averaging"))
				{
					if (options["radius_averaging"] == "arithmetic")
					{
						use_geometric_mean = false;
					}
					else if (options["radius_averaging"] == "geometric")
					{
						use_geometric_mean = true;
					}
					else
					{
						Log.warn() << "AmberNonBonded: unknown method for averaging LJ radii: '" 
								<< options["radius_averaging"] << "'. Using arithmetic mean." << std::endl;
					}
				}
			} 
			else if (hasVariable("A") && hasVariable("B"))
			{
				format_ = A_B_FORMAT;
			}		
			else if (hasVariable("alpha") && hasVariable("N") && hasVariable("R"))
			{
				// ?????
				format_ = SLATER_KIRKWOOD_FORMAT;
				Log.error() << "LennardJones::extractSection: Slater Kirkwood format not yet supported!" << endl;

				return false;
			}
							 
		}

		// build a two dimensional array of the atom types
		// loop variable
		Size	i;
		const AtomTypes&	atom_types = parameters.getAtomTypes();
		number_of_atom_types_ = atom_types.getNumberOfTypes();
		
		// allocate two onedimensional fields for the two parameters
		A_.resize(number_of_atom_types_);
		B_.resize(number_of_atom_types_);
		Aij_.resize(number_of_atom_types_ * number_of_atom_types_);
		Bij_.resize(number_of_atom_types_ * number_of_atom_types_);
		is_defined_.resize(number_of_atom_types_);

		for (i = 0; i < number_of_atom_types_; i++) 
		{
			is_defined_[i] = false;
		}

		// the indices of the columns containing the values
		Size index_A = 0;
		Size index_B = 0;
		if (format_ == A_B_FORMAT)
		{
			index_A = getColumnIndex("A");
			index_B = getColumnIndex("B");
		} 
		else if (format_ == EPSILON_R_FORMAT) 
		{
			index_A = getColumnIndex("epsilon");
			index_B = getColumnIndex("R");
		}
		else if (format_ == SLATER_KIRKWOOD_FORMAT) 
		{
			index_A = getColumnIndex("epsilon");
			index_B = getColumnIndex("R");
		}

		// try to identify the units of A and B
		// and set the two conversion factors
		double factor_A = 1.0;
		double factor_B = 1.0;
		
		if (format_ == A_B_FORMAT)
		{
			if (options.has("unit_A"))
			{ //?????
			}
		}
		// in EPSILON_R_FORMAT epsilon is in A and R is in B
		if (format_ == EPSILON_R_FORMAT)
		{
			if (options.has("unit_epsilon"))
			{
				if (options["unit_epsilon"] == "kcal/mol")
				{
					factor_A = Constants::JOULE_PER_CAL;
				}
				if (options["unit_epsilon"] == "cal/mol")
				{
					factor_A = Constants::JOULE_PER_CAL * 0.001;					
				}
				if (options["unit_epsilon"] == "J/mol")
				{
					factor_A = 0.001;
				}
			}

			if (options.has("unit_R"))
			{
				if (options["unit_R"] == "pm")
				{
					factor_B = 0.1;
				}
			}
		}

		Atom::Type		atom_type;
		String				key;
		for (i = 0; i < getNumberOfKeys(); ++i)
		{
			// get the key
			key = getKey(i);
			if (atom_types.hasType(key))
			{
				// get the two parameters
				atom_type = atom_types.getType(key);
				
				// retrieve the two values
				float A = getValue(i, index_A).toFloat() * factor_A;
				float B = getValue(i, index_B).toFloat() * factor_B;


				// store the values
				is_defined_[atom_type] = true;
				A_[atom_type] = A;
				B_[atom_type] = B;


				// check for the sign of the parameters: they have to be positive!
				if ((A < 0) || (B < 0))
				{
					if (format_ == EPSILON_R_FORMAT)
					{
						Log.warn() << "LennardJones::extractSection: VdW parameter may not be negative: type = " << atom_type << " (" << key << "), eps = " << A 
											 << ", r = " << B << endl;
					} 
					else 
					{
						Log.warn() << "LennardJones::extractSection: VdW parameter may not be negative: type = " << atom_type << " (" << key << "), A = " << A 
											 << ", B = " << B << endl;
					}
				}

			} 
			else	
			{
				Log.warn() << "LennardJones::extractSection: unknown atom type in Lennard Jones parameters: " << key << "   i = " << i << endl;
			}
		}

		// now assemble all Lennard Jones parameter for all known atom types
			
		for (i = 0; i < number_of_atom_types_; i++)
		{
			for (Size j = i; j < number_of_atom_types_; j++)
			{
				// calculate the two indices for the Aij/Bij fields
				Index index = (Index)(i + number_of_atom_types_ * j);
				Index sym_index = (Index)(j + number_of_atom_types_ * i);

				if (is_defined_[j] && is_defined_[i])
				{
					// calculate the values for A and B if in eps/R format
					if (format_ == EPSILON_R_FORMAT)
					{
						double R;
						if (!use_geometric_mean)
						{
							R = B_[i] + B_[j];
						}
						else
						{
							R = 2.0 * sqrt(B_[i] * B_[j]);
						}
						double R3 = R * R * R;
						double R6 = R3 * R3;
						double epsilon = sqrt(A_[i] * A_[j]);
						Aij_[index] = epsilon * R6 * R6;
						Bij_[index] = 2.0 * epsilon * R6;
					} 
					else 
					{
						// compute and assign Aij/Bij:
						// Aij = Ai * Aj, Aji = Aij
						// Bij = Bi * Bj, Bji = Bij
						Aij_[index] = A_[i] * A_[j];
						Bij_[index] = B_[i] * B_[j];
					}
					Aij_[sym_index] = Aij_[index];
					Bij_[sym_index] = Bij_[index];
				} 
				else	
				{
					Aij_[index] = 0.0;
					Bij_[index] = 0.0;
					Aij_[sym_index] = 0.0;
					Bij_[sym_index] = 0.0;
				}				
			}
		}

		return true;
	}


	bool LennardJones::hasParameters(Atom::Type I, Atom::Type J) const
		throw()
	{
		if (I < 0 || I >= (Index)number_of_atom_types_)
		{
			return false;
		}

		if (J < 0 || J >= (Index)number_of_atom_types_)
		{
			return false;
		}

		return (is_defined_[I] && is_defined_[J]);
	}


	LennardJones::Values LennardJones::getParameters(Atom::Type I, Atom::Type J) const throw()
	{
		LennardJones::Values parameters;
		assignParameters(parameters, I, J);
		return parameters;
	}


	bool LennardJones::assignParameters(LennardJones::Values& parameters, Atom::Type I, Atom::Type J) const 
		throw()
	{
		if (hasParameters(I, J)) 
		{
			parameters.A = Aij_[I * number_of_atom_types_ + J];
			parameters.B = Bij_[I * number_of_atom_types_ + J];
			
			return true;
		}

		return false;
	}


	bool LennardJones::operator == (const LennardJones& lj) const throw()
	{
		return (ParameterSection::operator == (lj)
						&& (A_ 		== lj.A_) 
						&& (B_ 		== lj.B_)
						&& (Aij_ 	== lj.Aij_) 
						&& (Bij_ 	== lj.Bij_));
	}
	 
} // namespace BALL
