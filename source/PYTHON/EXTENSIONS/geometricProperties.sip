// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// $Id: geometricProperties.sip,v 1.7 2004/02/25 10:47:12 oliver Exp $
//

class BoundingBoxProcessor
//  : public UnaryProcessor<Atom>
	: AtomProcessor
{
%TypeHeaderCode
	#include <BALL/STRUCTURE/geometricProperties.h>
	#include "sipBALLAtom.h"
	using namespace BALL;
	typedef UnaryProcessor<Atom> AtomProcessor;
%End
  public:
  virtual bool start() throw();
  virtual bool finish() throw();
  Processor::Result operator()(Atom& /atom/) throw();
%MethodCode
  sipRes = sipCpp->operator () (*a0);
%End
	SimpleBox3 getBox() const;
  const Vector3& getLower() const;
  const Vector3& getUpper() const;
};

class GeometricCenterProcessor
//	:	public UnaryProcessor<Atom> 
	: AtomProcessor
{
%TypeHeaderCode
	#include <BALL/STRUCTURE/geometricProperties.h>
	#include "sipBALLAtom.h"
	using namespace BALL;
	typedef UnaryProcessor<Atom> AtomProcessor;
%End
  public:
  virtual bool start() throw();
  virtual bool finish() throw();
  Processor::Result operator()(Atom& /atom/) throw();
%MethodCode
  sipRes = sipCpp->operator () (*a0);
%End
  Vector3& getCenter();
};

class FragmentDistanceCollector
//  : public UnaryProcessor<Composite> 
	: CompositeProcessor
{
%TypeHeaderCode
	#include <BALL/STRUCTURE/geometricProperties.h>
	#include "sipBALLAtom.h"
	using namespace BALL;
	typedef UnaryProcessor<Composite> CompositeProcessor;
%End
  public:
  FragmentDistanceCollector();
  FragmentDistanceCollector(const Composite& /composite/);
  FragmentDistanceCollector(const Composite& /composite/, float /distance/);
	~FragmentDistanceCollector() throw();
  virtual bool start() throw();
  virtual bool finish() throw();
  Processor::Result operator()(Composite& /composite/) throw();
%MethodCode
  sipRes = sipCpp->operator () (*a0);
%End
  Size getNumberOfFragments();
  void setComposite(const Composite& /composite/);
  const Composite* getComposite() const;
  float getDistance() const;
  void setDistance(float /distance/);
  // vector<Fragment*> fragments;
};


%ModuleHeaderCode
	#include <BALL/STRUCTURE/geometricProperties.h>
	using namespace BALL;
%End
Angle calculateTorsionAngle(const Atom& /a1/, const Atom& /a2/, const Atom& /a3/, const Atom& /a4/);
Angle calculateBondAngle(const Atom& /a1/, const Atom& /a2/, const Atom& /a3/); 
