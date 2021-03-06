// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// $Id: pierottiCavFreeEnergyProcessor.C,v 1.10 2002/02/27 12:24:06 sturm Exp $

#include <BALL/SOLVATION/pierottiCavFreeEnergyProcessor.h>
#include <BALL/STRUCTURE/numericalSAS.h>

using namespace std;

namespace BALL
{

	const char* PierottiCavFreeEnergyProcessor::Option::VERBOSITY = "verbosity";
	const char* PierottiCavFreeEnergyProcessor::Option::SOLVENT_NUMBER_DENSITY 
		= "solvent_number_density";
	const char* PierottiCavFreeEnergyProcessor::Option::PRESSURE = "pressure";
	const char* PierottiCavFreeEnergyProcessor::Option::ABSOLUTE_TEMPERATURE 
		= "absolute_temperature";
	const char* PierottiCavFreeEnergyProcessor::Option::PROBE_RADIUS 
		= "probe_radius";
	
	const int PierottiCavFreeEnergyProcessor::Default::VERBOSITY = 0;
	const float PierottiCavFreeEnergyProcessor::Default::SOLVENT_NUMBER_DENSITY 
		= 3.33253e-2;
	const float PierottiCavFreeEnergyProcessor::Default::PRESSURE = 1.01325e5;
	const float PierottiCavFreeEnergyProcessor::Default::ABSOLUTE_TEMPERATURE 
		= 298.0;
	const float PierottiCavFreeEnergyProcessor::Default::PROBE_RADIUS = 1.385;


	PierottiCavFreeEnergyProcessor::PierottiCavFreeEnergyProcessor() throw()
	{
		setDefaultOptions();

		valid_ = true;
	}


	PierottiCavFreeEnergyProcessor::PierottiCavFreeEnergyProcessor
	(const PierottiCavFreeEnergyProcessor& proc) throw()
		:	EnergyProcessor(proc)
	{
	}


	PierottiCavFreeEnergyProcessor::~PierottiCavFreeEnergyProcessor() throw()
	{
		clear();

		valid_ = false;
	}


	void PierottiCavFreeEnergyProcessor::clear() throw()
	{
		setDefaultOptions();

		valid_ = true;
	}

        const PierottiCavFreeEnergyProcessor& PierottiCavFreeEnergyProcessor::operator = (const PierottiCavFreeEnergyProcessor& proc) throw()     
        {
	         valid_=proc.valid_;
                 energy_=proc.energy_;
                 fragment_=proc.fragment_;  
                 return *this;
        }

        bool PierottiCavFreeEnergyProcessor::operator == (const PierottiCavFreeEnergyProcessor& proc) const throw()
        {
          bool result;
		if ((fragment_ == 0) && (proc.fragment_ == 0))
		{
			result = ((energy_ == proc.energy_) && (valid_ == proc.valid_));
		}
		else
		{
			if ((fragment_ == 0) || (proc.fragment_ == 0))
			{
				result = false;
			}
			else
			{
				result = ((*fragment_ == *proc.fragment_) 
						&& (energy_ 	 == proc.energy_)
						&& (valid_ 	 == proc.valid_));
			}
		}
		return result;
	}

	bool PierottiCavFreeEnergyProcessor::finish() throw()
	{

		// first check for user settings

		int verbosity = (int) options.getInteger(Option::VERBOSITY);
		// rho is the number density of the solvent (i. e. water) [1/m^3]
		double rho = options.getReal(Option::SOLVENT_NUMBER_DENSITY) * 1e30;
		// the pressure [ Pa ]
		double P = options.getReal(Option::PRESSURE);
		// the temperature
		double T = options.getReal(Option::ABSOLUTE_TEMPERATURE);
		// the solvent radius [ A ]
		double solvent_radius = options.getReal(Option::PROBE_RADIUS);
		if (verbosity > 0)
		{
			Log.info() << "Using a probe radius of " << solvent_radius << " A" <<
				endl;
		}
		
		// now compute some constant terms (names as in Pierotti, Chem. Rev.
		// 76(6):717--726, 1976)

		double s1 = 2 * solvent_radius * 1e-10; // [ m ]
		double s1_3 = s1 * s1 * s1; // [ m^3 ]
		double y = Constants::PI * s1_3 * (rho/6);	// [ 1 ]
		double y_frac = y/(1-y); // 
		double NkT = Constants::AVOGADRO * Constants::BOLTZMANN * T; // [ J/mol ]
		
		if (verbosity > 0)
		{
			Log.info() << "y = " << y << endl;
			Log.info() << "y_frac = " << y_frac << endl;
		}

		HashMap<const Atom*,float> atom_areas;
		calculateSASAtomAreas(*fragment_, atom_areas, solvent_radius);
		
		// R is the relation between solute radius and solvent radius [ 1 ]
		double R; 
		// deltaGspher is the cavitatonal energy of a spherical solute [ J/mol ]
		double deltaGspher; 
		// deltaGcav is the cavitatonal energy of the molecule [ J/mol ]
		double deltaGcav = 0;

		HashMap<const Atom*,float>::Iterator it = atom_areas.begin();
		for (; +it; ++it)
		{
			R = 2 * it->first->getRadius() * 1e-10 / s1;

			/* the following code implements the Pierotti '76 variant, considering
			 * relations of radii */

			deltaGspher =	-log(1.0 - y) 
					+ 3.0 * y_frac * R
					+ ( 3.0 * y_frac + 4.5 * y_frac * y_frac ) * ( R * R )
					+ y * P / ( rho * NkT ) * ( R * R * R );
			deltaGspher *= NkT;

			R = it->first->getRadius() * 1e-10 + s1 / 2.0;
			deltaGcav += it->second * 1e-20 /
				( 4 * Constants::PI * R * R ) * deltaGspher;
		}

		// return energy in units of kJ/mol
		energy_ = deltaGcav/1000;
		return 1;
	}


	void PierottiCavFreeEnergyProcessor::setDefaultOptions() throw()
	{
		options.setDefaultInteger(Option::VERBOSITY, Default::VERBOSITY);
		options.setDefaultReal(Option::SOLVENT_NUMBER_DENSITY, 
				Default::SOLVENT_NUMBER_DENSITY);
		options.setDefaultReal(Option::ABSOLUTE_TEMPERATURE,
				Default::ABSOLUTE_TEMPERATURE);
		options.setDefaultReal(Option::PRESSURE, Default::PRESSURE);
		options.setDefaultReal(Option::PROBE_RADIUS, Default::PROBE_RADIUS);
	}

} // namespace BALL
