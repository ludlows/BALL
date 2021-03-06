// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// $Id: amber.C,v 1.30 2005/01/18 21:46:49 amoll Exp $
//
// Molecular Mechanics: Amber force field class
//

#include <BALL/SYSTEM/path.h>
#include <BALL/MOLMEC/AMBER/amber.h>
#include <BALL/MOLMEC/AMBER/amberStretch.h>
#include <BALL/MOLMEC/AMBER/amberBend.h>
#include <BALL/MOLMEC/AMBER/amberTorsion.h>
#include <BALL/MOLMEC/AMBER/amberNonBonded.h>
#include <BALL/MOLMEC/COMMON/assignTypes.h>
#include <BALL/MOLMEC/PARAMETER/templates.h>

using namespace std;

namespace BALL 
{
	const char* AmberFF::Option::FILENAME = "filename";
	const char* AmberFF::Option::NONBONDED_CUTOFF = "nonbonded_cutoff";
	const char* AmberFF::Option::VDW_CUTOFF = "vdw_cutoff";
	const char* AmberFF::Option::VDW_CUTON = "vdw_cuton";
	const char* AmberFF::Option::ELECTROSTATIC_CUTOFF = "electrostatic_cutoff";
	const char* AmberFF::Option::ELECTROSTATIC_CUTON = "electrostatic_cuton";
	const char* AmberFF::Option::SCALING_VDW_1_4 = "SCAB";
	const char* AmberFF::Option::SCALING_ELECTROSTATIC_1_4 = "SCEE";
	const char* AmberFF::Option::DISTANCE_DEPENDENT_DIELECTRIC = "DDDC"; 
	const char* AmberFF::Option::ASSIGN_CHARGES = "assign_charges"; 
	const char* AmberFF::Option::ASSIGN_TYPENAMES = "assign_type_names"; 
	const char* AmberFF::Option::ASSIGN_TYPES = "assign_types"; 
	const char* AmberFF::Option::OVERWRITE_CHARGES = "overwrite_non-zero_charges"; 
	const char* AmberFF::Option::OVERWRITE_TYPENAMES = "overwrite_non-empty_typenames"; 

	const char* AmberFF::Default::FILENAME = "Amber/amber96.ini";
	const float AmberFF::Default::NONBONDED_CUTOFF = 20.0;
	const float AmberFF::Default::VDW_CUTOFF = 15.0;
	const float AmberFF::Default::VDW_CUTON = 13.0;
	const float AmberFF::Default::ELECTROSTATIC_CUTOFF = 15.0;
	const float AmberFF::Default::ELECTROSTATIC_CUTON = 13.0;
	const float AmberFF::Default::SCALING_ELECTROSTATIC_1_4 = 2.0;
	const float AmberFF::Default::SCALING_VDW_1_4 = 2.0;
	const bool  AmberFF::Default::DISTANCE_DEPENDENT_DIELECTRIC = false;   
	const bool	AmberFF::Default::ASSIGN_CHARGES = true;
	const bool	AmberFF::Default::ASSIGN_TYPENAMES = true;
	const bool	AmberFF::Default::ASSIGN_TYPES = true;
	const bool	AmberFF::Default::OVERWRITE_CHARGES = true;
	const bool	AmberFF::Default::OVERWRITE_TYPENAMES = false;

	// Default constructor
	AmberFF::AmberFF() 
		: ForceField(),
			filename_(Default::FILENAME),
			parameters_initialized_(false)
	{
		// set the force field name
		setName("Amber [" + filename_ + "]");

		// create the component list
		insertComponent(new AmberStretch(*this));
		insertComponent(new AmberBend(*this));
		insertComponent(new AmberTorsion(*this));
		insertComponent(new AmberNonBonded(*this));
	}

  // Constructor initialized with a system
  AmberFF::AmberFF(System& system)
    : ForceField(),
			filename_(Default::FILENAME),
			parameters_initialized_(false)
  {
		// create the component list
		insertComponent(new AmberStretch(*this));
		insertComponent(new AmberBend(*this));
		insertComponent(new AmberTorsion(*this));
		insertComponent(new AmberNonBonded(*this));

    bool result = setup(system);

		// set the force field name
		setName("Amber [" + filename_ + "]");

    if (!result)
    {
      Log.error() << " Force Field setup failed! " << endl;
      valid_ = false;
		}
	}

  // Constructor intialized with a system and a set of options
  AmberFF::AmberFF(System& system, const Options& new_options)
    : ForceField(),
			filename_(Default::FILENAME),
			parameters_initialized_(false)
  {
		// create the component list
		insertComponent(new AmberStretch(*this));
		insertComponent(new AmberBend(*this));
		insertComponent(new AmberTorsion(*this));
		insertComponent(new AmberNonBonded(*this));

    bool result = setup(system, new_options);

		// set the force field name (this has to be done after(!) setup,
		// otherwise filename_ is not yet set if it is used in options)
		setName("Amber [" + filename_ + "]");

    if (!result)
    {
      Log.error() << " Force Field setup failed! " << endl;
      valid_ = false;
		}
	}
 
	// copy constructor  
	AmberFF::AmberFF(const AmberFF& force_field)
		:	ForceField(force_field),
			filename_(force_field.filename_),
			parameters_initialized_(false)
	{
	}

	// destructor 
	AmberFF::~AmberFF()
	{
	}

	void AmberFF::clear()
		throw()
	{
		ForceField::clear();
		filename_ = Default::FILENAME;
		parameters_initialized_ = false;
	}

	const AmberFF& AmberFF::operator = (const AmberFF& force_field)
		throw()
	{
		// avoid self assignment
		if (&force_field != this)
		{
			ForceField::operator = (force_field);
			filename_ = force_field.filename_;
			parameters_initialized_ = force_field.parameters_initialized_;
		}
		
		return *this;
	}

	bool AmberFF::specificSetup()
		throw(Exception::TooManyErrors)
	{
		// check whether the system is assigned
		if (getSystem() == 0)
		{
			return false;
		}
 
		// check whether the parameter file name
		// is set in the options
		if (options.has(Option::FILENAME))
		{
			filename_ = options[Option::FILENAME];
			setName("Amber [" + filename_ + "]");
		} 
		else 
		{
			options[Option::FILENAME] = filename_;
		}

		// open parameter file
		Path    path;
		String  filename(path.find(filename_));

		if (filename == "") 
		{
			throw Exception::FileNotFound(__FILE__, __LINE__, filename_);
		}

		// initialize the force field parameters
		// and retrieve the atom types
		if (parameters_.getFilename() != filename || !parameters_initialized_)
		{
			parameters_.setFilename(filename);
			parameters_.init();

			// this is the first time, parameters was initialized
			// tell all components about it
			parameters_initialized_ = false;

			// retrieve global force field options
			ParameterSection global_options;
			global_options.extractSection(parameters_, "Options");
			// ?????: Iterator ueber global_options.options
			if (global_options.options.has(Option::NONBONDED_CUTOFF))
			{
				options.setDefault(Option::NONBONDED_CUTOFF, global_options.options[Option::NONBONDED_CUTOFF]);
			}
			if (global_options.options.has(Option::VDW_CUTOFF))
			{
				options.setDefault(Option::VDW_CUTOFF, global_options.options[Option::VDW_CUTOFF]);
			}
			if (global_options.options.has(Option::ELECTROSTATIC_CUTOFF))
			{
				options.setDefault(Option::ELECTROSTATIC_CUTOFF, global_options.options[Option::ELECTROSTATIC_CUTOFF]);
			}
			if (global_options.options.has(Option::SCALING_VDW_1_4))
			{
				options.setDefault(Option::SCALING_VDW_1_4, global_options.options[Option::SCALING_VDW_1_4]);
			}
			if (global_options.options.has(Option::SCALING_ELECTROSTATIC_1_4))
			{
				options.setDefault(Option::SCALING_ELECTROSTATIC_1_4, global_options.options[Option::SCALING_ELECTROSTATIC_1_4]);
			}
			if (global_options.options.has(Option::DISTANCE_DEPENDENT_DIELECTRIC))
			{
				options.setDefault(Option::DISTANCE_DEPENDENT_DIELECTRIC, global_options.options[Option::DISTANCE_DEPENDENT_DIELECTRIC]);
			}
		}

		// check the options whether types, type names, or charges 
		// have to be (re)assigned
		options.setDefaultBool(Option::ASSIGN_CHARGES, Default::ASSIGN_CHARGES);
		bool assign_charges = options.getBool(Option::ASSIGN_CHARGES);
		options.setDefaultBool(Option::ASSIGN_TYPES, Default::ASSIGN_TYPES);
		bool assign_types = options.getBool(Option::ASSIGN_TYPES);
		options.setDefaultBool(Option::ASSIGN_TYPENAMES, Default::ASSIGN_TYPENAMES);
		bool assign_type_names = options.getBool(Option::ASSIGN_TYPENAMES);
		options.setDefaultBool(Option::OVERWRITE_TYPENAMES, Default::OVERWRITE_TYPENAMES);
		bool overwrite_type_names = options.getBool(Option::OVERWRITE_TYPENAMES);
		options.setDefaultBool(Option::OVERWRITE_CHARGES, Default::OVERWRITE_CHARGES);
		bool overwrite_charges = options.getBool(Option::OVERWRITE_CHARGES);
		
		// extract template section (containing charges and atom types)
		if (assign_charges || assign_type_names)
		{
			Templates templates;
			templates.setMaximumUnassignedAtoms(max_number_of_errors_ - number_of_errors_);
			templates.extractSection(parameters_, "ChargesAndTypeNames");
			if (assign_charges && assign_type_names)
			{
				templates.assign(*getSystem(), overwrite_type_names, overwrite_charges);
			} 
			else 
			{
				if (assign_type_names)
				{
					templates.assignTypeNames(*getSystem(), overwrite_type_names);
				} 
				else 
				{
					templates.assignCharges(*getSystem(), overwrite_charges);
				}
			}
			
			HashSet<const Atom*>::ConstIterator it = templates.getUnassignedAtoms().begin();
			for (; it != templates.getUnassignedAtoms().end(); it++)
			{
			  getUnassignedAtoms().insert(*it);
			}

			number_of_errors_ += templates.getUnassignedAtoms().size();
	
			if (number_of_errors_ > max_number_of_errors_)
			{
				throw(Exception::TooManyErrors(__FILE__, __LINE__));
			}
		}
		if (assign_types)
		{
			// convert the type names to types
			AssignTypeProcessor type_proc(parameters_.getAtomTypes());
			type_proc.setMaximumUnassignedAtoms(max_number_of_errors_ - number_of_errors_);
			getSystem()->apply(type_proc);			

			HashSet<const Atom*>::ConstIterator it = type_proc.getUnassignedAtoms().begin();
			for (; it != type_proc.getUnassignedAtoms().end(); it++)
			{
			  getUnassignedAtoms().insert(*it);
			}

			number_of_errors_ += type_proc.getUnassignedAtoms().size();

			if (number_of_errors_ > max_number_of_errors_)
			{
				throw(Exception::TooManyErrors(__FILE__, __LINE__));
			}
		}

		return true;
	}

	Size AmberFF::getUpdateFrequency() const
	{
		return 20;
	}

	double AmberFF::getStretchEnergy() const
	{
		ForceFieldComponent* component = getComponent("Amber Stretch");
		if (component != 0)
		{
			return component->getEnergy();
		} 
		else 
		{
			return 0.0;
		}
	}

	double AmberFF::getBendEnergy() const
	{
		ForceFieldComponent* component = getComponent("Amber Bend");
		if (component != 0)
		{
			return component->getEnergy();
		} 
		else 
		{
			return 0;
		}
	}

	double AmberFF::getTorsionEnergy() const
	{
		ForceFieldComponent* component = getComponent("Amber Torsion");
		if (component != 0)
		{
			return component->getEnergy();
		} 
		else 
		{
			return 0;
		}
	}

	double AmberFF::getVdWEnergy() const
	{
		const ForceFieldComponent* component = getComponent("Amber NonBonded");
		if (component != 0)
		{
			const AmberNonBonded* nonbonded_component = dynamic_cast<const AmberNonBonded*>(component);
			if (nonbonded_component != 0)
			{
				return nonbonded_component->getVdwEnergy();
			}
		}

		return 0;
	}

	double AmberFF::getESEnergy() const
	{
		const ForceFieldComponent* component = getComponent("Amber NonBonded");
		if (component != 0)
		{
			const AmberNonBonded* nonbonded_component = dynamic_cast<const AmberNonBonded*>(component);
			if (nonbonded_component != 0)
			{
				return nonbonded_component->getElectrostaticEnergy();
			}
		}

		return 0;
	}

	double AmberFF::getNonbondedEnergy() const
	{
		const ForceFieldComponent* component = getComponent("Amber NonBonded");
		if (component != 0)
		{
			return component->getEnergy();
		}

		return 0;
	}

	bool AmberFF::hasInitializedParameters() const
	{
		return parameters_initialized_;
	}

	String AmberFF::getResults() const
		throw()
	{
		String result = String("\n")
		+ "AMBER Energy:\n"
		+ " - electrostatic     : " +String(getESEnergy())+  " kJ/mol\n" 
		+ " - van der Waals     : " +String(getVdWEnergy())+  " kJ/mol\n"
		+ " - bond stretch      : " +String(getStretchEnergy())+  " kJ/mol\n"
		+ " - angle bend        : " +String(getBendEnergy())+  " kJ/mol\n" 
		+ " - torsion           : " +String(getTorsionEnergy())+  " kJ/mol\n" 
		+ "---------------------------------------\n" 
		+ "  total energy       : " +String(getEnergy()) + " kJ/mol\n";
		return result;
	}
	
} // namespace BALL
