// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// $Id: atomBondModelBaseProcessor.sip,v 1.1.6.1 2005/11/11 17:21:18 leonhardt Exp $
//

class AtomBondModelBaseProcessor 
  : ModelProcessor
{
%TypeHeaderCode
	#include <BALL/VIEW/MODELS/atomBondModelBaseProcessor.h>
	using namespace BALL;
	using namespace BALL::VIEW;
%End
  
	public:
  
		AtomBondModelBaseProcessor() throw();
  	AtomBondModelBaseProcessor(const AtomBondModelBaseProcessor&) throw();
  	~AtomBondModelBaseProcessor() throw();
  	virtual void clear() throw();
  	void set(const AtomBondModelBaseProcessor&) throw();
		//  const AtomBondModelBaseProcessor& operator = (const AtomBondModelBaseProcessor& processor) throw();
  	//Processor::Result operator() (Composite&);
		//  void dump(std::ostream& s = std::cout, Size depth = 0) const throw();
		virtual void clearComposites() throw();
		virtual bool createGeometricObjects() throw();
};
