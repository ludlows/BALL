// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// $Id: surface.sip,v 1.7 2003/10/29 22:54:47 oliver Exp $
//

class Surface
{
%TypeHeaderCode
	#include <BALL/MATHS/surface.h>
	using namespace BALL;
%End
  public:
  class Triangle
 	{
		public:
	  Index v1;
	  Index v2;
	  Index v3;
  };

	public:
  Surface();
  Surface(const Surface& /surface/);
  ~Surface();
  void clear();
  void set(const Surface& /box/);
	//  const Surface& operator = (const Surface& box);
  void get(Surface& /box/) const;
  void readMSMSFile(const String& /vert_filename/, const String& /face_filename/);
  float getArea() const;
  bool operator == (const Surface& /box/) const;
  bool operator != (const Surface& /box/) const;
  // vector<Vector3>  vertex;
  // vector<Vector3>  normal;
  // vector<Triangle> triangle;

	Size getNumberOfTriangles() const throw();
	Size getNumberOfVertices() const throw();
	Size getNumberOfNormals() const throw();
	Triangle& getTriangle(Position /index/) throw();
	const Triangle& getTriangle(Position /index/) const throw();
	Vector3& getVertex(Position /index/) throw();
	const Vector3& getVertex(Position /index/) const throw();
	Vector3& getNormal(Position /index/) throw();
	const Vector3& getNormal(Position /index/) const throw();
	void clearVertices();
	void clearTriangles();
	void clearNormals();
	void pushBackTriangle(const Triangle& /t/) throw();
	void pushBackVertex(const Vector3& /v/) throw();
	void pushBackNormal(const Vector3& /n/) throw();
	void resizeTriangles(Size /size/);
	void resizeNormals(Size /size/);
	void resizeVertices(Size /size/);
};
