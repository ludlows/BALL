// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// $Id: Bruker2DFile_test.C,v 1.2 2003/06/02 14:34:29 oliver Exp $
//

#include <BALL/CONCEPT/classTest.h>

///////////////////////////

#include <BALL/FORMAT/bruker2DFile.h>

///////////////////////////

START_TEST(Bruker2DFile, "$Id: Bruker2DFile_test.C,v 1.2 2003/06/02 14:34:29 oliver Exp $")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace BALL;

Bruker2DFile* ptr = 0;
CHECK(Bruker2DFile::Bruker2DFile() throw())
	ptr = new Bruker2DFile;
	TEST_NOT_EQUAL(ptr, 0)
RESULT											

CHECK(Bruker2DFile::~Bruker2DFile() throw())
	delete ptr;
RESULT											

/// ????

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
