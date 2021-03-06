// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// $Id: PeakList_test.C,v 1.8 2003/11/10 08:47:45 anhi Exp $
#include <BALL/CONCEPT/classTest.h>

///////////////////////////

#include <BALL/NMR/peakList.h>

///////////////////////////

START_TEST(PeakList, "$Id: PeakList_test.C,v 1.8 2003/11/10 08:47:45 anhi Exp $")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace BALL;


/// PeakList1D
PeakList1D* pl_ptr = 0;
CHECK(PeakList1D::PeakList())
	pl_ptr = new PeakList1D;
	TEST_NOT_EQUAL(pl_ptr, 0)
RESULT


CHECK(PeakList1D::~PeakList())
	delete pl_ptr;
RESULT


PeakList1D pl;
PeakList1D::PeakType p1;
PeakList1D::PeakType p2;
p1.setPosition(1.0);
p2.setPosition(4.0);
p1.setIntensity(2.0);
p2.setIntensity(8.0);
pl.push_back(p1);
pl.push_back(p2);

CHECK(PeakList1D::PeakList(const PeakList& peak_list))
	TEST_EQUAL(pl.size(), 2)
	PeakList1D pl_copy(pl);
	TEST_EQUAL(pl_copy.size(), 2)
	TEST_REAL_EQUAL(pl_copy.front().getPosition(), 1.0)
	TEST_REAL_EQUAL(pl_copy.back().getPosition(), 4.0)
	TEST_REAL_EQUAL(pl_copy.front().getIntensity(), 2.0)
	TEST_REAL_EQUAL(pl_copy.back().getIntensity(), 8.0)
RESULT


CHECK(PeakList1D::scale(float x))
	PeakList1D pl2(pl);
	pl2.scale(2.0);
	TEST_EQUAL(pl2.size(), 2)
	TEST_REAL_EQUAL(pl2.front().getPosition(), 1.0)
	TEST_REAL_EQUAL(pl2.back().getPosition(), 4.0)
	TEST_REAL_EQUAL(pl2.front().getIntensity(), 4.0)
	TEST_REAL_EQUAL(pl2.back().getIntensity(), 16.0)
RESULT


CHECK(PeakList1D::getMaxIntensity() const )
	PeakList1D pl2;
	TEST_REAL_EQUAL(pl2.getMaxIntensity(), -Limits<float>::max())
	TEST_REAL_EQUAL(pl.getMaxIntensity(), 8.0)
RESULT


CHECK(PeakList1D::getMinIntensity() const )
	PeakList1D pl2;
	TEST_REAL_EQUAL(pl2.getMinIntensity(), Limits<float>::max())
	TEST_REAL_EQUAL(pl.getMinIntensity(), 2.0)
RESULT


/// PeakList2D
PeakList2D* pl2D_ptr = 0;
CHECK(PeakList2D::PeakList())
	pl2D_ptr = new PeakList2D;
	TEST_NOT_EQUAL(pl2D_ptr, 0)
RESULT


CHECK(PeakList2D::~PeakList())
	delete pl2D_ptr;
RESULT


PeakList2D pl_2D;
PeakList2D::PeakType p_2D_1;
PeakList2D::PeakType p_2D_2;
p_2D_1.setPosition(Vector2(1.0, 2.0));
p_2D_2.setPosition(Vector2(4.0, 5.0));
p_2D_1.setIntensity(2.0);
p_2D_2.setIntensity(8.0);
pl_2D.push_back(p_2D_1);
pl_2D.push_back(p_2D_2);

CHECK(PeakList2D::PeakList(const PeakList& peak_list))
	TEST_EQUAL(pl_2D.size(), 2)
	PeakList2D pl_2D_copy(pl_2D);
	TEST_EQUAL(pl_2D_copy.size(), 2)
	TEST_REAL_EQUAL(pl_2D_copy.front().getPosition().x, 1.0)
	TEST_REAL_EQUAL(pl_2D_copy.front().getPosition().y, 2.0)
	TEST_REAL_EQUAL(pl_2D_copy.back().getPosition().x, 4.0)
	TEST_REAL_EQUAL(pl_2D_copy.back().getPosition().y, 5.0)
	TEST_REAL_EQUAL(pl_2D_copy.front().getIntensity(), 2.0)
	TEST_REAL_EQUAL(pl_2D_copy.back().getIntensity(), 8.0)
RESULT


CHECK(PeakList2D::scale(float x))
	PeakList2D pl_2D_2(pl_2D);
	pl_2D_2.scale(2.0);
	TEST_EQUAL(pl_2D_2.size(), 2)
	TEST_REAL_EQUAL(pl_2D_2.front().getPosition().x, 1.0)
	TEST_REAL_EQUAL(pl_2D_2.front().getPosition().y, 2.0)
	TEST_REAL_EQUAL(pl_2D_2.back().getPosition().x, 4.0)
	TEST_REAL_EQUAL(pl_2D_2.back().getPosition().y, 5.0)
	TEST_REAL_EQUAL(pl_2D_2.front().getIntensity(), 4.0)
	TEST_REAL_EQUAL(pl_2D_2.back().getIntensity(), 16.0)
RESULT


CHECK(PeakList2D::getMaxIntensity() const )
	PeakList2D pl_2D_2;
	TEST_REAL_EQUAL(pl_2D_2.getMaxIntensity(), -Limits<float>::max())
	TEST_REAL_EQUAL(pl_2D.getMaxIntensity(), 8.0)
RESULT


CHECK(PeakList2D::getMinIntensity() const )
	PeakList2D pl_2D_2;
	TEST_REAL_EQUAL(pl_2D_2.getMinIntensity(), Limits<float>::max())
	TEST_REAL_EQUAL(pl_2D.getMinIntensity(), 2.0)
RESULT


/// PeakList3D
PeakList3D* pl3D_ptr = 0;
CHECK(PeakList3D::PeakList())
	pl3D_ptr = new PeakList3D;
	TEST_NOT_EQUAL(pl3D_ptr, 0)
RESULT


CHECK(PeakList3D::~PeakList())
	delete pl3D_ptr;
RESULT


PeakList3D pl_3D;
PeakList3D::PeakType p_3D_1;
PeakList3D::PeakType p_3D_2;
p_3D_1.setPosition(Vector3(1.0, 2.0, 3.0));
p_3D_2.setPosition(Vector3(4.0, 5.0, 6.0));
p_3D_1.setIntensity(2.0);
p_3D_2.setIntensity(8.0);
pl_3D.push_back(p_3D_1);
pl_3D.push_back(p_3D_2);

CHECK(PeakList3D::PeakList(const PeakList& peak_list))
	TEST_EQUAL(pl_3D.size(), 2)
	PeakList3D pl_3D_copy(pl_3D);
	TEST_EQUAL(pl_3D_copy.size(), 2)
	TEST_REAL_EQUAL(pl_3D_copy.front().getPosition().x, 1.0)
	TEST_REAL_EQUAL(pl_3D_copy.front().getPosition().y, 2.0)
	TEST_REAL_EQUAL(pl_3D_copy.front().getPosition().z, 3.0)
	TEST_REAL_EQUAL(pl_3D_copy.back().getPosition().x, 4.0)
	TEST_REAL_EQUAL(pl_3D_copy.back().getPosition().y, 5.0)
	TEST_REAL_EQUAL(pl_3D_copy.back().getPosition().z, 6.0)
	TEST_REAL_EQUAL(pl_3D_copy.front().getIntensity(), 2.0)
	TEST_REAL_EQUAL(pl_3D_copy.back().getIntensity(), 8.0)
RESULT


CHECK(PeakList3D::scale(float x))
	PeakList3D pl_3D_2(pl_3D);
	pl_3D_2.scale(2.0);
	TEST_EQUAL(pl_3D_2.size(), 2)
	TEST_REAL_EQUAL(pl_3D_2.front().getPosition().x, 1.0)
	TEST_REAL_EQUAL(pl_3D_2.front().getPosition().y, 2.0)
	TEST_REAL_EQUAL(pl_3D_2.front().getPosition().z, 3.0)
	TEST_REAL_EQUAL(pl_3D_2.back().getPosition().x, 4.0)
	TEST_REAL_EQUAL(pl_3D_2.back().getPosition().y, 5.0)
	TEST_REAL_EQUAL(pl_3D_2.back().getPosition().z, 6.0)
	TEST_REAL_EQUAL(pl_3D_2.front().getIntensity(), 4.0)
	TEST_REAL_EQUAL(pl_3D_2.back().getIntensity(), 16.0)
RESULT


CHECK(PeakList3D::getMaxIntensity() const )
	PeakList3D pl_3D_2;
	TEST_REAL_EQUAL(pl_3D_2.getMaxIntensity(), -Limits<float>::max())
	TEST_REAL_EQUAL(pl_3D.getMaxIntensity(), 8.0)
RESULT


CHECK(PeakList3D::getMinIntensity() const )
	PeakList3D pl_3D_2;
	TEST_REAL_EQUAL(pl_3D_2.getMinIntensity(), Limits<float>::max())
	TEST_REAL_EQUAL(pl_3D.getMinIntensity(), 2.0)
RESULT


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

END_TEST

