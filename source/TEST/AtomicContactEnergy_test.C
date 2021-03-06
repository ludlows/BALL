// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// $Id: AtomicContactEnergy_test.C,v 1.9 2003/05/07 18:23:39 oliver Exp $
//

#include <BALL/CONCEPT/classTest.h>

///////////////////////////

#include <BALL/ENERGY/atomicContactEnergy.h>
#include <BALL/FORMAT/PDBFile.h>

///////////////////////////

START_TEST(AtomicContactEnergy, "$Id: AtomicContactEnergy_test.C,v 1.9 2003/05/07 18:23:39 oliver Exp $")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace BALL;

CHECK(calculateACE())
	// read protein A (trypsin from 2ptc)
	STATUS("reading PDB file data/ACE_test_A.pdb")
	PDBFile pdb_file_A("data/ACE_test_A.pdb");
	TEST_EQUAL(pdb_file_A.good(), true)
	System A;
	pdb_file_A >> A;
	pdb_file_A.close();
	TEST_EQUAL(A.countAtoms(), 1629)

	// read protein B (BPTI from 2ptc)
	STATUS("reading PDB file data/ACE_test_B.pdb")
	PDBFile pdb_file_B("data/ACE_test_B.pdb");
	TEST_EQUAL(pdb_file_B.good(), true)
	System B;
	pdb_file_B >> B;
	pdb_file_B.close();
	TEST_EQUAL(B.countAtoms(), 454)

	// calculate the atomic contact energy of A and B
	PRECISION(1e-2)
	STATUS("calculating the ACE of A and B")
	float ACE_A = calculateACE(A);
	TEST_REAL_EQUAL(ACE_A, -2222.15)
	float ACE_B = calculateACE(B);
	TEST_REAL_EQUAL(ACE_B, -562.36)

	// join A and B and calculate the total energy
	STATUS("calculating the ACE of AB")
	A.splice(B);
	float ACE_AB = calculateACE(A);
	TEST_REAL_EQUAL(ACE_AB, -2806.2)

	// calculate the difference in binding energies
	float dG_bind = ACE_AB - ACE_A - ACE_B;
	TEST_REAL_EQUAL(dG_bind, -21.6891)
RESULT											

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
