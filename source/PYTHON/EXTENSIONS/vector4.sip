// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// $Id: vector4.sip,v 1.8 2004/09/22 21:07:06 oliver Exp $
//

class Vector4
{
%TypeHeaderCode
	#include <BALL/MATHS/vector4.h>
	#include "sipBALLIndex.h"
  using namespace BALL;
  using namespace BALL::Exception;
%End
  public:
  // BALL_CREATE_NODEEP(Vector4)
  Vector4();
  // Vector4(float* /ptr/);
  Vector4(float /x/, float /y/, float /z/, float /h/ = 1);
  Vector4(const Vector4& /vector/);
  ~Vector4() 
		throw();
  // void set(float* /ptr/);
  void set(float /rx/, float /ry/, float /rz/, float /rh/ = 1);
  void set(const Vector4& /vector/);
  // Vector4& operator = (float* ptr);
  // Vector4& operator = (const Vector4& vector);
  // void get(float* /ptr/) const;
  void get(float /rx/, float /ry/, float /rz/, float /rh/) const;
  void get(Vector4& /vector/) const;
  void swap(Vector4& /vector/);
  float getLength() const;
  float getSquareLength() const;
  Vector4& normalize();
  static const Vector4& getZero();
  static const Vector4& getUnit();
  void set(float /value/ = 1);
	float operator [] (Index /index/);
%MethodCode
	sipRes = sipCpp->operator [] (*a0);
%End

  // float operator [] (Index /index/) const;
  // ???? Vector4 operator + () const;
  // ???? Vector4 operator - () const;
  Vector4& operator += (const Vector4& /vector/);
  Vector4& operator -= (const Vector4& /vector/);
  Vector4 operator * (float /scalar/);
  Vector4& operator *= (float /scalar/);
  Vector4 operator / (float /scalar/);
  Vector4& operator /= (float /scalar/);
  float operator * (const Vector4& /vector/) const;
  float getDistance(const Vector4& /vector/) const;
  float getSquareDistance(const Vector4& /vector/) const;
  bool operator == (const Vector4& /vector/) const;
  bool operator != (const Vector4& /vector/) const;
  bool isOrthogonalTo(Vector4& /vector/) const;
  bool isValid() const;
  // void dump(std::ostream& s = std::cout, Size depth = 0) const;

  float x;
  float y;
  float z;
  float h;
};
 
