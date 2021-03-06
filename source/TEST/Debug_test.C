// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// $Id: Debug_test.C,v 1.1 2004/02/24 08:19:42 oliver Exp $
//

#include <BALL/CONFIG/config.h>
// Make sure the debug flag is set, independent of what 
// configure says
#define BALL_DEBUG
#include <BALL/COMMON/debug.h>

#include <BALL/CONCEPT/classTest.h>

///////////////////////////

START_TEST(COMMON_debug, "$Id: Debug_test.C,v 1.1 2004/02/24 08:19:42 oliver Exp $")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace BALL;

CHECK(BALL_PRECONDITION_EXCEPTION(condition, message))
	TEST_EXCEPTION(Exception::Precondition, BALL_PRECONDITION_EXCEPTION(1 == 0, "test1"))
	BALL_PRECONDITION_EXCEPTION(1 == 1, "test2")
RESULT

CHECK(BALL_POSTCONDITION_EXCEPTION(condition, message))
	TEST_EXCEPTION(Exception::Postcondition, BALL_POSTCONDITION_EXCEPTION(1 == 0, "test1"))
	BALL_POSTCONDITION_EXCEPTION(1 == 1, "test2")
RESULT

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
