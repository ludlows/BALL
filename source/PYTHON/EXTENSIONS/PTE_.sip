// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// $Id: PTE_.sip,v 1.5 2003/11/01 09:27:36 oliver Exp $
//

class Element
  : PropertyManager
{
%TypeHeaderCode
	#include <BALL/KERNEL/PTE.h>
	using namespace BALL;
%End
   public:
   Element();
   Element(const Element& /element/);
//   Element(const Element&, bool = true); this is nonsense
//   Element
//    (const String& name,
//     const String& symbol,
//     Group group,
//     Period period,
//     AtomicNumber atomic_umber,
//     float atomic_weight,
//     float atomic_radius,
//     float covalent_radius,
//     float van_der_waals_radius,
//     float electronegativity);
//   virtual ~Element();
//   static const char *getErrorMessage(ErrorCode );
   void setName(const String&);
   const String& getName() const;
   void setSymbol(const String&);
   const String& getSymbol() const;
   void setGroup(int);
   int getGroup() const;
   void setPeriod(int);
   int getPeriod() const;
   void setAtomicNumber(int);
   int getAtomicNumber() const;
   void setAtomicWeight(float );
   float getAtomicWeight() const;
   void setAtomicRadius(float);
   float getAtomicRadius() const;
   void setCovalentRadius(float);
   float getCovalentRadius() const;
   void setVanDerWaalsRadius(float);
   float getVanDerWaalsRadius() const;
   void setElectronegativity(float);
   float getElectronegativity() const;
   bool operator == (const Element&) const;
	 bool operator !=(const Element&) const;
   bool operator <(const Element&) const;
   bool operator <=(const Element&) const;
   bool operator >=(const Element&) const;
   bool operator >(const Element&) const;
   //friend ostream& operator << (ostream& s, const Element&);
   bool isUnknown() const;

	 static Element UNKNOWN;
 };

class PTE_
	: PropertyManager
{
	public:
	PTE_();
	PTE_(const PTE_&);
	//   PTE_(const PTE_&, bool);
	//   virtual ~PTE_();
	static Element& getElement(int);
	static Element& getElement(const String&);
	// const Element& operator [] (const String&);

	const Element& operator [] (Index /index/) const;
%MethodCode
  sipRes = &sipCpp->operator [] (*a0);
%End

	// static bool apply(UnaryProcessor<Element>& applicator);
};
