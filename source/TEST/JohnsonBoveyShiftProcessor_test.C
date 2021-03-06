// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// $Id: JohnsonBoveyShiftProcessor_test.C,v 1.13 2004/11/18 19:34:15 oliver Exp $
//

#include <BALL/CONCEPT/classTest.h>

///////////////////////////

#include <BALL/NMR/johnsonBoveyShiftProcessor.h>
#include <BALL/FORMAT/HINFile.h>
#include <BALL/FORMAT/PDBFile.h>

///////////////////////////

START_TEST(JohnsonBoveyShiftProcessor, "$Id: JohnsonBoveyShiftProcessor_test.C,v 1.13 2004/11/18 19:34:15 oliver Exp $")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace BALL;
using namespace std;

JohnsonBoveyShiftProcessor* sp = 0;
CHECK(JohnsonBoveyShiftProcessor::JohnsonBoveyShiftProcessor() throw())
	sp = new JohnsonBoveyShiftProcessor;
	TEST_NOT_EQUAL(sp, 0)
RESULT


CHECK(JohnsonBoveyShiftProcessor::~JohnsonBoveyShiftProcessor() throw())
  delete sp;
RESULT

Parameters parameters("data/JohnsonBoveyShiftProcessor_test.ini");
HINFile f("data/JohnsonBoveyShiftProcessor_test.hin");
System S;
f >> S;
f.close();


CHECK(JohnsonBoveyShiftProcessor::JohnsonBoveyShiftProcessor(const JohnsonBoveyShiftProcessor& processor) throw())
	JohnsonBoveyShiftProcessor sp1;
	TEST_EQUAL(sp1.isValid(), false)
	sp1.setParameters(parameters);
	sp1.init();
	TEST_EQUAL(sp1.isValid(), true)
	JohnsonBoveyShiftProcessor sp2(sp1);
	TEST_EQUAL(sp2.isValid(), true)
RESULT


CHECK(JohnsonBoveyShiftProcessor::init() throw())
  JohnsonBoveyShiftProcessor sp;
	sp.setParameters(parameters);
	TEST_EQUAL(sp.isValid(), false)
	sp.init();
	TEST_EQUAL(sp.isValid(), true)
RESULT


CHECK(JohnsonBoveyShiftProcessor::start() throw())
	JohnsonBoveyShiftProcessor sp;
	TEST_EQUAL(sp.start(), false)
	sp.setParameters(parameters);
	TEST_EQUAL(sp.start(), false)
	sp.init();	
	TEST_EQUAL(sp.start(), true)
RESULT


CHECK(JohnsonBoveyShiftProcessor::finish() throw())
	// tested below
RESULT


CHECK(JohnsonBoveyShiftProcessor::Processor::Result operator () (Composite& composite) throw())
	// tested below
RESULT


CHECK(chemical shifts/without rings)
	PRECISION(0.0001)
	JohnsonBoveyShiftProcessor sp;
	sp.setParameters(parameters);
	sp.init();
	TEST_EQUAL(sp.isValid(), true)
	TEST_EQUAL(S.countAtoms(), 31)
	
	if (S.countAtoms() == 31)
	{
		S.apply(sp);

		AtomIterator atom_it = S.beginAtom();
		Position i = 0;
		for (; +atom_it; ++atom_it)
		{
			if (atom_it->hasProperty(JohnsonBoveyShiftProcessor::PROPERTY__RING_CURRENT_SHIFT))
			{
				i++;
				TEST_EQUAL(atom_it->getProperty(JohnsonBoveyShiftProcessor::PROPERTY__RING_CURRENT_SHIFT).getFloat(), 0.0)
			}
		}
		TEST_EQUAL(i, 15)
	}	
RESULT


S.destroy();
f.open("data/JohnsonBoveyShiftProcessor_test2.hin");
f >> S;
CHECK(chemical shifts/with rings)
	StringHashMap<float> rc_shifts;
	ifstream infile("data/JohnsonBoveyShiftProcessor_test.dat");
	String name;
	float shift;
	while (infile.good())
	{	
		infile >> name >> shift;
		if (name != "")
		{
			rc_shifts.insert(name, shift);
		}
	}
	TEST_EQUAL(rc_shifts.size(), 160)

	JohnsonBoveyShiftProcessor sp;
	sp.setParameters(parameters);
	sp.init();
	TEST_EQUAL(sp.isValid(), true)
	TEST_EQUAL(S.countAtoms(), 328)

	PRECISION(0.01)
	
	if (S.countAtoms() == 328)
	{
		S.apply(sp);

		AtomIterator atom_it = S.beginAtom();
		Position i = 0;
		for (; +atom_it; ++atom_it)
		{
			if (atom_it->hasProperty(JohnsonBoveyShiftProcessor::PROPERTY__RING_CURRENT_SHIFT))
			{
				STATUS("atom " << atom_it->getFullName() << " has shift property = " << atom_it->getProperty(JohnsonBoveyShiftProcessor::PROPERTY__RING_CURRENT_SHIFT).getFloat())
				shift = atom_it->getProperty(JohnsonBoveyShiftProcessor::PROPERTY__RING_CURRENT_SHIFT).getFloat();
				if (shift != 0.0)
				{
					STATUS("shift of " << atom_it->getFullName() << ": " << shift)
					TEST_EQUAL(rc_shifts.has(atom_it->getFullName()), true)
					if (rc_shifts.has(atom_it->getFullName()))
					{
						TEST_REAL_EQUAL(shift, rc_shifts[atom_it->getFullName()])
						i++;
					}
				}
			}
		}
		TEST_EQUAL(i, 160)
	}	
RESULT


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
