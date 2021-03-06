// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// $Id: RingPerceptionProcessor_test.C,v 1.2 2004/08/07 19:38:25 oliver Exp $
//

#include <BALL/CONCEPT/classTest.h>

///////////////////////////

#include <BALL/QSAR/ringPerceptionProcessor.h>
#include <BALL/FORMAT/SDFile.h>
#include <BALL/KERNEL/system.h>
#include <BALL/KERNEL/atom.h>
#include <BALL/KERNEL/bond.h>
#include <BALL/KERNEL/forEach.h>
#include <BALL/KERNEL/molecule.h>

///////////////////////////
START_TEST(RingPerceptionProcessor, "$Id: RingPerceptionProcessor_test.C,v 1.2 2004/08/07 19:38:25 oliver Exp $")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace BALL;

SDFile infile("data/descriptors_test.sdf");
System S;
infile >> S;
infile.close();

Molecule * molecule;
Size limit = S.countMolecules();

CHECK(RingPerceptionProcessor)	

	RingPerceptionProcessor rpp;
	S.apply(rpp);

	Size results_atoms[9] = {6, 6, 6, 6, 6, 13, 0, 6, 4};
	Size results_bonds[9] = {6, 6, 6, 6, 6, 14, 0, 6, 6};

	for (Size i = 0; i != limit; ++i)
	{
		molecule = S.getMolecule(i);
		Size num_ringed = 0;

		// atoms
		AtomConstIterator a_it = molecule->beginAtom();
		for (; a_it != molecule->endAtom(); ++a_it)
		{
			if (a_it->getProperty("InRing").getBool())
			{		
				++num_ringed;
			}
		}
		TEST_EQUAL(num_ringed, results_atoms[i])

		// bonds
		num_ringed = 0;
		a_it = molecule->beginAtom();
		Atom::BondConstIterator b_it = a_it->beginBond();
		BALL_FOREACH_BOND(*molecule, a_it, b_it)
		{
			if (b_it->getProperty("InRing").getBool())
			{
				++num_ringed;
			}
		}
		TEST_EQUAL(num_ringed, results_bonds[i])
	}
		
RESULT

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
