// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// $Id: HBondModel.sip,v 1.1.6.1 2005/11/11 17:21:18 leonhardt Exp $
//

class HBondModelProcessor  
	: AtomBondModelBaseProcessor
{
%TypeHeaderCode
	#include <BALL/VIEW/MODELS/HBondModel.h>
	using namespace BALL;
	using namespace BALL::VIEW;
%End
  public:
  
	HBondModelProcessor() throw();
  HBondModelProcessor(const HBondModelProcessor&) throw();
  ~HBondModelProcessor() throw();
  virtual void clear() throw();
  void set(const HBondModelProcessor&) throw();
//  const HBondModelProcessor& operator = (const HBondModelProcessor& model) throw();
  //Processor::Result operator() (Composite&);
  void setRadius(float) throw();
  float getRadius() const throw();
};
