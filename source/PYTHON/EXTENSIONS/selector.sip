// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// $Id: selector.sip,v 1.5 2003/11/01 09:27:40 oliver Exp $
//
class Selector
//  : UnaryProcessor<Composite>
	: CompositeProcessor
{
%TypeHeaderCode
	#include <BALL/KERNEL/selector.h>
	using namespace BALL;
	typedef UnaryProcessor<Composite> CompositeProcessor;
%End
  public:
  // BALL_CREATE_NODEEP(Selector)
  Selector();
  Selector(const Selector& /selector/);
  Selector(const String& /expression_string/);
  ~Selector() throw();
  Processor::Result operator () (Composite& /composite/) throw();
%MethodCode
	sipRes = sipCpp->operator () (*a0);
%End

  virtual bool start() throw();
  Size getNumberOfSelectedAtoms() const;
};
