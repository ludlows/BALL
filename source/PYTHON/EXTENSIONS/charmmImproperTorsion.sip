// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// $Id: charmmImproperTorsion.sip,v 1.5 2004/12/27 17:24:56 amoll Exp $
//

class CharmmImproperTorsion 
  : ForceFieldComponent
{
%TypeHeaderCode
	#include <BALL/MOLMEC/CHARMM/charmmImproperTorsion.h>
	#include <BALL/COMMON/exception.h>
	using namespace BALL;
	using namespace BALL::Exception;
%End
  public:
  // BALL_CREATE(CharmmImproperTorsion)
  CharmmImproperTorsion();
  CharmmImproperTorsion(ForceField& /force_field/);
  CharmmImproperTorsion(const CharmmImproperTorsion& /charmm_stretch/);
  ~CharmmImproperTorsion();
  virtual bool setup() throw(TooManyErrors);
  virtual double updateEnergy();
  virtual void updateForces();
};
  
