// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// $Id: atomContainer.sip,v 1.7 2003/11/01 09:27:36 oliver Exp $
//

class AtomContainer
  : Composite, 
		PropertyManager
{
%TypeHeaderCode
	#include <BALL/KERNEL/atomContainer.h>
	using namespace BALL;
%End

  public:
  enum Property
  {
    NUMBER_OF_PROPERTIES = 0
  };
  AtomContainer();
  AtomContainer(const AtomContainer&, bool = true);
  AtomContainer(const String&);
	~AtomContainer() throw();
  virtual void clear() throw();
  virtual void destroy() throw();
	//  virtual void persistentWrite(PersistenceManager& pm, const char* = 0) const;
	//  virtual void persistentRead(PersistenceManager& pm);
  void set(const AtomContainer&, bool = true);
	//  AtomContainer& operator = (const AtomContainer&);
  void get(AtomContainer&, bool = true) const;
  void swap(AtomContainer&);
  void setName(const String&);
  const String& getName() const;
  AtomContainer* getSuperAtomContainer();
  const AtomContainer* getSuperAtomContainer() const;
  AtomContainer* getAtomContainer(Position);
  const AtomContainer* getAtomContainer(Position) const;
  Atom* getAtom(Position);
  const Atom* getAtom(Position) const;
  Atom* getAtom(const String&);
  const Atom* getAtom(const String&) const;
  Size countAtomContainers() const;
  Size countAtoms() const;
  Size countBonds() const;
  Size countInterBonds() const;
  Size countIntraBonds() const;
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
  void destroyBonds();
  bool isSubAtomContainerOf(const AtomContainer&) const;
  bool isSuperAtomContainerOf(const AtomContainer&) const;
  virtual bool isValid() const throw();
	//  virtual void dump(ostream& s = cout,unsigned long depth = 0) const;
	//  virtual void read(istream& s);
	//  virtual void write(ostream& s) const;

	// workaround for template method apply()...
	bool apply(AtomProcessor& /processor/);
	bool apply(CompositeProcessor& /processor/);
	bool apply(FragmentProcessor& /processor/);
	bool apply(ResidueProcessor& /processor/);

	
  SIP_PYOBJECT __str__();
%MethodCode
	sipRes = PyString_FromString(String(String("AtomContainer ") + sipCpp->getName() 
						+ " { " + String(sipCpp->countAtoms()) + " atoms }").c_str());
%End
};
