// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// $Id: gradient.sip,v 1.3 2003/10/29 22:54:41 oliver Exp $
//

class Gradient
{
%TypeHeaderCode
#include <BALL/MOLMEC/COMMON/gradient.h>
using namespace BALL;
%End
 public:

  Gradient();
  Gradient(const AtomVector&);
  Gradient(const Gradient&, bool);

	// Gradient& operator = (const Gradient& rhs);
  //Gradient& operator = (const AtomVector& rhs);
  void set(const Gradient&);
  void set(const AtomVector&);
  void negate();
  // double operator * (const Gradient& gradient) const;
  void normalize();
  int size() const ;
  // using vector<Vector3>::operator [];
  void invalidate();
  bool isValid() const;
  // using vector<Vector3>::begin;
  // using vector<Vector3>::end;
  double norm;
  double inv_norm;
  double rms;
};
