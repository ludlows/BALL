// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// $Id: charmmTorsion.sip,v 1.5 2004/12/27 17:24:57 amoll Exp $
//

class CharmmTorsion 
  : ForceFieldComponent
{
%TypeHeaderCode
	#include <BALL/MOLMEC/CHARMM/charmmTorsion.h>
	#include <BALL/COMMON/exception.h>
	using namespace BALL;
	using namespace BALL::Exception;
%End

  public:
//  struct SingleCharmmTorsion 
//  {
//   Atom* atom1;
//   Atom* atom2;
//   Atom* atom3;
//   Atom* atom4;
//   float  V;
//   unsigned char f;
//   float  phase;
//   SingleCharmmTorsion()
//    : atom1(0),
//     atom2(0),
//     atom3(0),
//     atom4(0),
//     V(0),
//     f(0),
//     phase(0)
//   {
//   }
//   SingleCharmmTorsion(CosineTorsion::SingleData& t)
//   {
//    atom1 = t.atom1;
//    atom2 = t.atom2;
//    atom3 = t.atom3;
//    atom4 = t.atom4;
//    V = t.values.V / t.values.n;
//    f = (unsigned char)t.values.f;
//    phase = (BALL::Constants::PI / 180.0) * t.values.phase;
//   }
//  };
  // BALL_CREATE(CharmmTorsion)
  CharmmTorsion();
  CharmmTorsion(ForceField& /force_field/);
  CharmmTorsion(const CharmmTorsion& /charmm_stretch/);
  ~CharmmTorsion();
  virtual bool setup() throw(TooManyErrors);
  virtual double updateEnergy();
  virtual void updateForces();
};
 
