// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// $Id: geometricTransformations.sip,v 1.6 2004/02/25 10:47:13 oliver Exp $
//

class TranslationProcessor
//  : public UnaryProcessor<Atom> 
	: AtomProcessor
{
%TypeHeaderCode
	#include <BALL/STRUCTURE/geometricTransformations.h>
	#include "sipBALLAtom.h"
	using namespace BALL;
	typedef UnaryProcessor<Atom> AtomProcessor;
%End
  public:
  TranslationProcessor();
  TranslationProcessor(const Vector3& /translation/);
  void setTranslation(const Vector3& /translation/);
  const Vector3& getTranslation() const;
  Processor::Result operator()(Atom& /atom/);
%MethodCode
  sipRes = sipCpp->operator () (*a0);
%End
};

class TransformationProcessor
//  :public UnaryProcessor<Atom> 
	: AtomProcessor
{
%TypeHeaderCode
	#include <BALL/STRUCTURE/geometricTransformations.h>
	#include "sipBALLAtom.h"
	using namespace BALL;
	typedef UnaryProcessor<Atom> AtomProcessor;
%End
  public:
  TransformationProcessor();
  TransformationProcessor(const Matrix4x4& /transformation/);
  void setTransformation(const Matrix4x4& /transformation/);
  const Matrix4x4& getTransformation() const;
  Processor::Result operator()(Atom& /atom/);
%MethodCode
  sipRes = sipCpp->operator () (*a0);
%End
};
