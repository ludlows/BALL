// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// $Id: box3.sip,v 1.8 2003/11/01 09:27:37 oliver Exp $
//

class Box3
{
%TypeHeaderCode
	#include <BALL/MATHS/box3.h>
	using namespace BALL;
%End
	public:
  Box3();
  Box3(Box3&);
  Box3(Vector3&, Vector3&);

  void set(Box3&);
  // Box3& operator = (Box3& box);
  float getSurface() const;
  float getVolume() const;
  float getWidth() const;
  float getHeight() const;
  float getDepth() const;
  bool operator == (const Box3& /box/) const;
  bool operator != (Box3& /box/) const;
  bool isValid() const;
  // void dump(std::ostream& s = std::cout, Size depth = 0) const;


  SIP_PYOBJECT __str__();
%MethodCode
	static String tmp;
	tmp = "Box3 { (";
	tmp += String(sipCpp->getWidth()) + " x ";
	tmp += String(sipCpp->getHeight()) + " x ";
	tmp += String(sipCpp->getDepth()) + ") @ ( ";
	tmp += String(sipCpp->getPoint().x) + " ";
	tmp += String(sipCpp->getPoint().y) + " ";
	tmp += String(sipCpp->getPoint().z) + ") }";
	sipRes = PyString_FromString(tmp.c_str());
%End
};
