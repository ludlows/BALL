// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// $Id: forceFieldComponent.sip,v 1.5 2004/12/27 17:24:57 amoll Exp $
//

class ForceFieldComponent 
{
%TypeHeaderCode
	#include <BALL/MOLMEC/COMMON/forceFieldComponent.h>
	#include <BALL/COMMON/exception.h>
	using namespace BALL;
	using namespace BALL::Exception;
%End
  public:
  ForceFieldComponent();
  ForceFieldComponent(ForceField& /force_field/);
  ForceFieldComponent(const ForceFieldComponent& /force_field/);
  virtual bool setup() throw(TooManyErrors);
  ForceField* getForceField() const;
  void setForceField(ForceField& /force_field/);
  void setName(const String&);
  String getName() const;
  virtual double getEnergy() const;
  virtual double updateEnergy();
  virtual void updateForces();
};
