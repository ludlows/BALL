// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// $Id: quaternion.sip,v 1.6 2003/10/29 22:54:45 oliver Exp $
//

class Quaternion
{
%TypeHeaderCode
	#include <BALL/MATHS/quaternion.h>
	using namespace BALL;
%End
  public:
  // BALL_CREATE_NODEEP(Quaternion)
  Quaternion();
  Quaternion(const Quaternion& /q/);
  Quaternion(const Vector3& /axis/, float /angle/);
  Quaternion(float /x/, float /y/, float /z/, float /angle/);
  ~Quaternion();
  void set(const Quaternion& /q/);
  void set(const Vector3& /axis/, float /angle/);
  void set(float /x/, float /y/, float /z/, float /angle/);
  // Quaternion& operator = (const Quaternion& q);
  void get(Quaternion& /q/) const;
  void setIdentity();
  void swap(Quaternion& /q/);
  float getAngle() const;
  Vector3 getAxis();
  Matrix4x4& getRotationMatrix(Matrix4x4& /m/) const;
  // ???? Quaternion operator - () const;
  Quaternion getInverse() const;
  Quaternion getConjugate() const;
  Quaternion& operator += (const Quaternion& /q/);
  Quaternion& operator -= (const Quaternion& /q/);
  bool operator == (const Quaternion& /q/) const;
  bool operator != (const Quaternion& /q/) const;
  // void dump(std::ostream& s = std::cout, Size depth = 0) const;
  float i;
  float j;
  float k;
  float angle;
};
