// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// $Id: periodicBoundary.sip,v 1.4 2003/10/31 23:15:49 oliver Exp $
//

class PeriodicBoundary
{
%TypeHeaderCode
#include <BALL/MOLMEC/COMMON/periodicBoundary.h>
using namespace BALL;
%End
 public:
  PeriodicBoundary();
  PeriodicBoundary(const ForceField&);
  PeriodicBoundary(const PeriodicBoundary&);

  // PeriodicBoundary& operator = (const PeriodicBoundary&);
  bool setup();
  int generateMoleculesVector();
  void enable();
  void disable();
  SimpleBox3 getBox() const;
  void setBox(const SimpleBox3&);
  int addSolvent(const String&) const;
  int removeSolvent() const;
  bool isEnabled() const;
  void updateMolecules();
  Options* options;
};
