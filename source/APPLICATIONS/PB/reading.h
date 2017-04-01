// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// $Id: reading.h,v 1.4 2002/02/27 12:20:32 sturm Exp $

enum RuleType
{
	CHARGES_AND_RADII,
	CHARGES,
	RADII
};

// read the contents of a PDB file and assign
// charges and radii 
void readPDBFile(const String& filename);

// read the contents of a PDB file and assign
// charges and radii 
void readHINFile(const String& filename);

// read the contents of a PDB file and assign radii 
void readHINFileNoAssignment(const String& filename);

// read an option file
void readOptionFile(const String& filename);

// read a radius file
void readRadiusFile(const String& filename);

// read a charge file
void readChargeFile(const String& filename);

// read a charge and radius rule file
void readRuleFile(const String& filename, RuleType rule_type);

// dump the positions, charges, and radii to a file
void dumpFile();
