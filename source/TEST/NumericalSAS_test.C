// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// $Id: NumericalSAS_test.C,v 1.8 2002/02/27 12:24:40 sturm Exp $
#include <BALL/CONCEPT/classTest.h>

///////////////////////////
#include <BALL/STRUCTURE/numericalSAS.h>
#include <BALL/KERNEL/fragment.h>
#include <BALL/MATHS/surface.h>
#include <BALL/DATATYPE/hashMap.h>
///////////////////////////

START_TEST(NumericalSAS, "$Id: NumericalSAS_test.C,v 1.8 2002/02/27 12:24:40 sturm Exp $")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace BALL;

CHECK(calculateSASArea(const AtomContainer&, float probe_radius, Size number_of_points))
	Fragment	f;
	Atom a1, a2;
	a1.setRadius(1.0);
	a2.setRadius(1.0);
	a2.setPosition(Vector3(10.0, 0.0, 0.0));

	f.insert(a1);
	f.insert(a2);

	float area = calculateSASArea(f, 1.5, 624);

	PRECISION(0.001)
	TEST_REAL_EQUAL(area, 157.07963)
RESULT

CHECK(calculateSASVolume(const AtomContainer&, float probe_radius, Size number_of_points))
	Fragment	f;
	Atom a1, a2;
	a1.setRadius(1.0);
	a2.setRadius(1.0);
	a2.setPosition(Vector3(10.0, 0.0, 0.0));

	f.insert(a1);
	f.insert(a2);

	float volume = calculateSASVolume(f, 1.5, 624);

	PRECISION(0.001)
	TEST_REAL_EQUAL(volume, 130.899)
RESULT

CHECK(calculateSASAtomAreas())
	Fragment	f;
	Atom a1, a2;
	a1.setRadius(1.0);
	a2.setRadius(1.0);
	a2.setPosition(Vector3(10.0, 0.0, 0.0));

	f.insert(a1);
	f.insert(a2);

	HashMap<const Atom*, float>	atom_map;

	float area = calculateSASAtomAreas(f, atom_map, 1.5, 624);

	PRECISION(0.001)
	TEST_REAL_EQUAL(area, 157.07963)
	TEST_REAL_EQUAL(atom_map[&a1], area / 2.0)
	TEST_REAL_EQUAL(atom_map[&a2], area / 2.0)
RESULT

CHECK(calculateSASPoints())
	Fragment	f;
	Atom a1, a2;
	a1.setRadius(1.0);
	a2.setRadius(1.0);
	a2.setPosition(Vector3(10.0, 0.0, 0.0));

	f.insert(a1);
	f.insert(a2);

	Surface surface;
	float area = calculateSASPoints(f, surface, 1.5, 624);

	PRECISION(0.001)
	TEST_REAL_EQUAL(area, 157.07963)
	TEST_EQUAL(surface.vertex.size(), 1284)
	TEST_EQUAL(surface.normal.size(), 1284)

	// sum up all normals to check for integrality of the 
	// surface elements
	float surface_elements = 0;
	for (Position i = 0; i < surface.normal.size(); i++)
	{
		surface_elements += surface.normal[i].getLength();
	}
	TEST_REAL_EQUAL(surface_elements, area)
RESULT

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
