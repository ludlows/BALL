// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// $Id: object.sip,v 1.6 2003/10/29 22:08:39 oliver Exp $
//

class Object
{
%TypeHeaderCode
#include <BALL/CONCEPT/object.h>
using namespace BALL;
%End
  public:
  Object() throw();
  Object(const Object&) throw();
  ~Object() throw();
  virtual void clear() throw();
  Size getHandle() const throw();
  static Size getNextHandle() throw();
  static Size getNewHandle() throw();
  bool operator == (const Object& /object/) const;
  bool operator != (const Object& /object/) const;
  bool operator < (const Object& /object/) const;
  bool operator <= (const Object& /object/) const;
 	bool operator >= (const Object& /object/) const;
  bool operator > (const Object& /object/) const;
  int compare(const Object&) const;
  // virtual void read(::std::istream& s);
  // virtual void write(::std::ostream& s) const;
  // friend ::std::istream& operator >> (::std::istream& s, Object& /object/);
  // friend ::std::ostream& operator << (::std::ostream& s, const Object& /object/);
  virtual bool isValid() const throw();
  //virtual void dump(::std::ostream& s = std::cout, unsigned long depth = 0) const;
};
