// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// $Id: label.sip,v 1.4 2004/05/27 19:49:59 oliver Exp $
//

class Label 
  : GeometricObject,
    Vertex
{
%TypeHeaderCode
	#include <BALL/VIEW/PRIMITIVES/label.h>
	using namespace BALL;
	using namespace BALL::VIEW;
%End
  public:
//  void dump(std::ostream& s = std::cout, Size depth = 0) const throw();
//  const Label& operator = (const Label& label) throw();
  Label() throw();
  Label(const Label&) throw();
  String getText() const throw();
  bool isValid() const throw();
  void clear() throw();
  void set(const Label&) throw();
  void setText(const String&) throw();
  void swap(Label&) throw();
  ~Label() throw();
};
