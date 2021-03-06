// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// $Id: RDFSection_test.C,v 1.6 2002/02/27 12:24:49 sturm Exp $
#include <BALL/CONCEPT/classTest.h>

///////////////////////////

#include <BALL/STRUCTURE/RDFSection.h>

///////////////////////////

START_TEST(RDFSection, "$Id: RDFSection_test.C,v 1.6 2002/02/27 12:24:49 sturm Exp $")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace BALL;

RDFSection* pointer;

CHECK(RDFSection::RDFSection())
  pointer = new RDFSection();
	TEST_NOT_EQUAL(pointer, 0)
RESULT


CHECK(RDFSection::~RDFSection())
  delete pointer;
RESULT


CHECK(RDFSection::RDFSection(const RDFSection& rdf_section))
  //?????
RESULT


CHECK(RDFSection::destroy())
  //?????
RESULT


CHECK(RDFSection::clear())
  //?????
RESULT


CHECK(RDFSection::set(const RDFSection& rdf_section))
  //?????
RESULT


CHECK(RDFSection::RDFSection& operator = (const RDFSection& rdf_section))
  //?????
RESULT


CHECK(RDFSection::extractSection(Parameters& parameters, const String& section_name))
  //?????
RESULT


CHECK(RDFSection::getRDF() const )
  //?????
RESULT




/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
