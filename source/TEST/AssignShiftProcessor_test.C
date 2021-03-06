// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// $Id: AssignShiftProcessor_test.C,v 1.16 2004/05/27 19:50:00 oliver Exp $
//

#include <BALL/CONCEPT/classTest.h>

#include <BALL/NMR/assignShiftProcessor.h>
#include <BALL/STRUCTURE/fragmentDB.h>
#include <BALL/STRUCTURE/residueChecker.h>
#include <BALL/FORMAT/HINFile.h>
#include <BALL/FORMAT/NMRStarFile.h>
#include <BALL/SYSTEM/path.h>

START_TEST(AssignShiftProcessor, "$Id: AssignShiftProcessor_test.C,v 1.16 2004/05/27 19:50:00 oliver Exp $")

using namespace BALL;
using namespace std;

System s;
CHECK(preparations)
	STATUS("opening HIN file")
	HINFile f;
	f.open("data/AssignShiftProcessor_test.hin");
	STATUS("reading structure")
	f >> s;
RESULT

FragmentDB db("fragments/Fragments.db");

CHECK(AssignShiftProcessor::AssignShiftProcessor(const vector<NMRAtomData*>& atom_data)/apply())
	NMRStarFile rs("data/AssignShiftProcessor_test.str");
	TEST_EQUAL(rs.getData()[0].atom_data.size(), 797)
	TEST_EQUAL(rs.getNumberOfAtoms(), 797)
	TEST_EQUAL(s.countAtoms(), 1944)

	AssignShiftProcessor asp(rs.getData()[0].atom_data);
	asp.setFragmentDB(&db);
	s.apply(asp);
RESULT

CHECK(results)
	Size number_of_shiftatoms = 0;
	AtomIterator atom_it = s.beginAtom();
	for (; +atom_it; ++atom_it)
	{
		//Log.error() << atom_it->getProperty(ShiftModule::PROPERTY__SHIFT).getFloat() << std::endl;
		if (atom_it->hasProperty(ShiftModule::PROPERTY__SHIFT))
		{
			number_of_shiftatoms++;
		}
	}
	TEST_EQUAL(number_of_shiftatoms, 940)
RESULT

END_TEST
