// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// $Id: molecule.sip,v 1.8 2003/11/15 10:59:06 oliver Exp $
//

class Molecule
  : AtomContainer
{
%TypeHeaderCode
	#include <BALL/KERNEL/molecule.h>
	using namespace BALL;
%End
  public:
  Molecule();
  Molecule(const Molecule&, bool = true);
  Molecule(const String&);
	~Molecule() throw();
  virtual void clear() throw();
  virtual void destroy() throw();
	//  void persistentWrite(PersistenceManager& pm, const char* name = 0) const;
	//  void persistentRead(PersistenceManager& pm);
  void set(const Molecule&, bool = true);
	//  Molecule& operator = (const Molecule&);
  void get(Molecule&, bool = true) const;
  void swap(Molecule&);
  System* getSystem();
  const System* getSystem() const;
  void prepend(Atom&);
  void append(Atom&);
  void insert(Atom&);
  void insertBefore(Atom&, Composite&);
  void insertAfter(Atom&, Composite&);
  bool remove(Atom&);
  void prepend(AtomContainer&);
  void append(AtomContainer&);
  void insert(AtomContainer&);
  void insertBefore(AtomContainer&, Composite&);
  void insertAfter(AtomContainer&, Composite&);
  void spliceBefore(AtomContainer&);
  void spliceAfter(AtomContainer&);
  void splice(AtomContainer&);
  bool remove(AtomContainer&);
  virtual bool isValid() const throw();
	//  virtual void dump(ostream& s = cout, unsigned long depth = 0) const;
	//  virtual void read(istream& s);
	//  virtual void write(ostream& s) const;

  SIP_PYOBJECT __repr__();
%MethodCode
	sipRes = PyString_FromString(String(String("Molecule ") + sipCpp->getName() 
						+ " { " + String(sipCpp->countAtoms()) + " atoms }").c_str());
%End

};
