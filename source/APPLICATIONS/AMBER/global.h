// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// $Id: global.h,v 1.3 2003/12/01 07:35:11 oliver Exp $
//

#include <BALL/STRUCTURE/fragmentDB.h>
#include <BALL/SOLVATION/poissonBoltzmann.h>
#include <BALL/STRUCTURE/defaultProcessors.h>
#include <BALL/MOLMEC/COMMON/radiusRuleProcessor.h>
#include <BALL/MOLMEC/COMMON/chargeRuleProcessor.h>

using namespace BALL;

// a pointer to the fragment data base (a pointer since we do
// not need it every time an want to save the loading time)
extern FragmentDB*           frag_db;

// the system that holds everything
extern System S;

// flags:

// the three main options:

// true, if a single point calculation is requested
extern bool energy_calculation;

// true, if an energy minimization is requested
extern bool energy_minimization;

// true, if the residue checker ist to be invoked
extern bool check_structures;

// true, if a selection string was given
extern bool use_selection;

//  true, if timing information and final options should be printed
extern bool verbose;

//  true, if the names should be normalized upon read
extern bool normalize_names;

//  true, if the bonds should be built upon read
extern bool build_bonds;

// the selection string
extern String selection;

// the optimizer gradient criterion
extern double max_gradient;

// the maximum number of iterations
extern Size max_iterations;

// the filename of the force field parameter file
extern String FF_filename;

// the filename of the output file
extern String out_filename;

//
extern bool sd_minimizer;

// The energy limit for warnings.
extern double energy_limit;
