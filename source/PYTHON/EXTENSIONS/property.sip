// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// $Id: property.sip,v 1.4 2003/10/29 22:54:43 oliver Exp $
//

class NamedProperty
  : PersistentObject
{
%TypeHeaderCode
	#include <BALL/CONCEPT/property.h>
	using namespace BALL;
%End
  public:
  NamedProperty(String&);
  NamedProperty(String&, bool);
  //NamedProperty(String&, int);
  //NamedProperty(String&, float);
  NamedProperty(String&, String&);
  NamedProperty(String&, PersistentObject&);
  NamedProperty(NamedProperty&);
  ~NamedProperty();
//    virtual void persistentWrite(PersistenceManager& pm, const char* name = "") const;
//    virtual void persistentRead(PersistenceManager& pm);
  int getType() const;
  String getName() const;
  bool getBool() const;
  int getInt() const;
  float getFloat() const;
  int getUnsignedInt() const;
  PersistentObject* getObject() const;
  String getString() const;
  NamedProperty();
};

class PropertyManager
{
  public:
  PropertyManager();
  PropertyManager(const PropertyManager&);
	~PropertyManager();
  void clear();
  void destroy();
  void set(const PropertyManager&);
	//   PropertyManager& operator = (const PropertyManager&);
  void get(PropertyManager&) const;
  void swap(PropertyManager&);
	//   BitVector& getBitVector();
	//   const BitVector& getBitVector() const;
	//   operator BitVector& ();
  void setProperty(int);
  void clearProperty(int);
  void toggleProperty(int);
  int countProperties() const;
  void setProperty(NamedProperty&);
  void setProperty(const String&);
  //void setProperty(const String&, bool);
  //void setProperty(const String&, int);
	//   void setProperty(const String&, unsigned int);
  void setProperty(const String&, float);
  void setProperty(const String&, String&);
  void setProperty(const String&, PersistentObject&);
  const NamedProperty& getProperty(const String&) const;
  void clearProperty(const String&);
	Size countNamedProperties() const;
  const NamedProperty& getNamedProperty(Position /index/) const
     throw(IndexOverflow);
  NamedProperty& getNamedProperty(Position /index/)
     throw(IndexOverflow);
  bool hasProperty(int) const;
  bool hasProperty(const String&) const;
	// friend ostream& operator << (ostream& s, const PropertyManager&);
	// friend istream& operator >> (istream& s, PropertyManager&);
	// void write(PersistenceManager& pm) const;
	// bool read(PersistenceManager& pm);
  bool isValid() const;
	// void dump(ostream& s = cout, unsigned long depth = 0) const;
};
