// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// $Id: common.sip,v 1.3 2003/10/29 22:54:38 oliver Exp $
//

%TypeHeaderCode
	#include <BALL/MATHS/common.h>
	using namespace BALL::Maths;
%End

// float abs(float /t/);
// float frac(float /t/);
// float max(float /a/, float /b/);
// float max(float /a/, float /b/, float /c/);
// float min(float /a/, float /b/);
// float min(float /a/, float /b/, float /c/);
float round(float /t/);
float sgn(float /t/);
bool isFinite(float /t/);
bool isInfinite(float /t/);
bool isNan(float /t/);
bool isZero(float /t/);
bool isNotZero(float /t/);
bool isEqual(float /a/, float /b/);
bool isNotEqual(float /a/, float /b/);
bool isLess(float /a/, float /b/);
bool isLessOrEqual(float /a/, float /b/);
bool isGreaterOrEqual(float /a/, float /b/);
bool isGreater(float /a/, float /b/);
long floor(float /t/);
long ceiling(float /t/);
Index compare(float /a/, float /b/);
bool isNear(float /a/, float /b/, float /max_diff/);
