// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// $Id: ruleProcessor.sip,v 1.5 2003/11/01 09:27:40 oliver Exp $
//

class RuleProcessor
  : AtomProcessor
{
%TypeHeaderCode
	#include <BALL/MOLMEC/COMMON/ruleProcessor.h>
%End
  public:
  // BALL_CREATE(RuleProcessor)
  RuleProcessor();
  RuleProcessor(INIFile& /file/, const String& /prefix/);
  RuleProcessor(const RuleProcessor& /rule_processor/);
  ~RuleProcessor();
  void clear();
  void destroy();   
  bool initialize(INIFile& /file/, const String& /prefix/);
  // const RuleProcessor& operator = (const RuleProcessor& /rule_processor/);
  void set(const RuleProcessor& /rule_processor/);
  virtual bool start();
  virtual bool finish();
  Processor::Result operator () (Atom& /atom/);
%MethodCode
  sipRes = sipCpp->operator () (*a0);
%End

  String evaluate(const Atom& /atom/);
  bool isValid() const;
  // void dump(std::ostream& s = std::cout) const;
};
