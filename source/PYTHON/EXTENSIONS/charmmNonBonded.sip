// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// $Id: charmmNonBonded.sip,v 1.5 2004/12/27 17:24:56 amoll Exp $
//

class CharmmNonBonded 
  : ForceFieldComponent
{
%TypeHeaderCode
	#include <BALL/MOLMEC/CHARMM/charmmNonBonded.h>
	#include <BALL/COMMON/exception.h>
	using namespace BALL;
	using namespace BALL::Exception;
%End
  public:
  // BALL_CREATE(CharmmNonBonded)
  CharmmNonBonded();
  CharmmNonBonded(ForceField& /force_field/);
  CharmmNonBonded(const CharmmNonBonded& /charmm_non_bonded/);
  ~CharmmNonBonded() throw();
  virtual bool setup() throw(TooManyErrors);
  virtual double updateEnergy() throw();
  virtual void updateForces() throw();
  virtual double getElectrostaticEnergy() const throw();
  virtual double getVdwEnergy() const throw();
  virtual double getSolvationEnergy() const throw();
};
