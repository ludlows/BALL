// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// $Id: global.h,v 1.9 2002/02/27 12:20:31 sturm Exp $

#include <BALL/STRUCTURE/fragmentDB.h>
#include <BALL/SOLVATION/poissonBoltzmann.h>
#include <BALL/STRUCTURE/defaultProcessors.h>
#include <BALL/MOLMEC/COMMON/radiusRuleProcessor.h>
#include <BALL/MOLMEC/COMMON/chargeRuleProcessor.h>

using namespace BALL;

// a pointer to the fragment data base (a pointer since we do
// not need it every time an want to save the loading time)
extern FragmentDB*           frag_db;

// the Finite Difference Poisson-Boltzmann object
extern FDPB                  fdpb;

// the options for the calculation
extern Options								options;

// charge and radius assignment processors
extern AssignChargeProcessor charges;
extern AssignRadiusProcessor radii;
extern ClearChargeProcessor	clear_charge_proc;

// rule-based assignment processors
extern ChargeRuleProcessor charge_rules;
extern RadiusRuleProcessor radius_rules;

// the system that holds everything
extern System S;

// the name of the dump file
extern String dump_file;

// the default probe sphere radius for surface calculations
extern float probe_radius;

// flags:

// the three main options:

// true, if the solvent excluded surface is to be calculated
extern bool ses_calculation;

// true, if the solvent accessible surface is to be calculated
extern bool sas_calculation;

// true, if a FDPB calculation is to ber performed
extern bool fdpb_calculation;

// the results of the SES and SAS area calculations
extern float total_SAS_area;
extern float total_SAS_volume;
extern float total_SES_area;
extern float total_SES_volume;

// a hash map containing the atom surfaces of the SAS (if calculated)
extern HashMap<const Atom*, float> surface_map;

//  true, if timing information and final options should be printed
extern bool verbose;

//  true, if the charges should be cleared upon read
extern bool clear_charges;

//  true, if the names should be normalized upon read
extern bool normalize_names;

//  true, if the bonds should be built upon read
extern bool build_bonds;

// true, if the rule-based assignment is to be used
extern bool use_charge_rules;
extern bool use_radius_rules;

// true, if the solvation free energy is to be calculated
// (difference calculation solvent/vacuum)
extern bool calculate_solvation_energy;



