// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
#include <BALL/MOLMEC/AMBER/amber.h>
#include <BALL/SOLVATION/poissonBoltzmann.h>
#include <BALL/FORMAT/HINFile.h>
#include <BALL/STRUCTURE/defaultProcessors.h>

using namespace BALL;
using namespace std;

int main()
{
	// read a system of molecules from a HyperChem file
	HINFile infile("test.hin");
	System s;
	infile >> s;
	
	// create an AMBER force field without non-bonded interactions
	AmberFF FF(s);
	FF.removeComponent("Amber NonBonded");
	
	// calculate the total energy
	float internal_energy = FF.updateEnergy();
	cout << "internal energy: " << internal_energy << " kJ/mol" << endl;

	
	// assign atom radii
	AssignRadiusProcessor radius_processor("/KM/comp-bio/BALL-data/radii/PARSE.siz");
	s.apply(radius_processor);

	// calculate the electrostatic part of the solvation energy
	FDPB fdpb(s);
	fdpb.solve();
	float solvent_energy = fdpb.getEnergy();
	
	fdpb.options[FDPB::Option::SOLVENT_DC] = 1.0;
	fdpb.setup(s);
	fdpb.solve();
	float vacuum_energy = fdpb.getEnergy();
	cout << "electrostatic solvation free energy: " 
		<< solvent_energy - vacuum_energy << endl;
	
	cout << "total energy: " << internal_energy - vacuum_energy + solvent_energy << " kJ/mol" << endl;
}
