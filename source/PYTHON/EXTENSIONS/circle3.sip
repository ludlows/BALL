// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// $Id: circle3.sip,v 1.4 2003/10/29 22:54:38 oliver Exp $
//

class Circle3
{
%TypeHeaderCode
	#include <BALL/MATHS/circle3.h>
	using namespace BALL;
%End
  public:
  // BALL_CREATE_NODEEP(Circle3)
  Circle3();
  Circle3(const Circle3& /circle/);
  Circle3(const Vector3& /point/, const Vector3& /normal/, float /radius/);
  ~Circle3();
  void swap(Circle3& /circle/);
  void set(const Circle3& /circle/);
  void set(const Vector3& /point/, const Vector3& /normal/, float /rad/);
  // Circle3& operator = (const Circle3& /circle/);
  void get(Circle3& /circle/) const;
  void get(Vector3& /point/, Vector3& /normal/, float /rhs/) const;
  bool operator == (const Circle3& /circle/) const;
  bool operator != (const Circle3& /circle/) const;
  bool has(const Vector3& /point/, bool /on_surface/ = false) const;
  bool isValid() const;
  // void dump(std::ostream& s = std::cout, Size depth = 0) const;
};
