// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// $Id: XYZFile_test.C,v 1.7 2003/07/03 13:30:41 oliver Exp $
//

#include <BALL/CONCEPT/classTest.h>

///////////////////////////

#include <BALL/FORMAT/XYZFile.h>
#include <BALL/KERNEL/system.h>
#include <BALL/KERNEL/PTE.h>

///////////////////////////

START_TEST(XYZFile, "$Id: XYZFile_test.C,v 1.7 2003/07/03 13:30:41 oliver Exp $")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace BALL;


XYZFile*	file_ptr;
CHECK(XYZFile::XYZFile())
	file_ptr = new XYZFile;
	TEST_NOT_EQUAL(file_ptr, 0)
RESULT


CHECK(XYZFile::~XYZFile())
	delete file_ptr;
RESULT


CHECK(XYZFile::XYZFile(const String& filename, File::OpenMode open_mode))
	// checking for default mode: reading
	XYZFile f("data/XYZFile_test1.xyz");
	TEST_EQUAL(f.isValid(), true)
RESULT


CHECK(XYZFile::read(System& system))
  XYZFile f("data/XYZFile_test1.xyz");
	System S;
	f.read(S);
	f.close();

	TEST_EQUAL(S.countAtoms(), 2)
	TEST_EQUAL(S.getName(), "TEST TEST")
	AtomIterator it = S.beginAtom();
	TEST_EQUAL(it->getElement().getSymbol(), "H")
	TEST_REAL_EQUAL(it->getPosition().x, 0.0)
	TEST_REAL_EQUAL(it->getPosition().y, 1.0)
	TEST_REAL_EQUAL(it->getPosition().z, 2.0)
	it++;
	TEST_EQUAL(it->getElement().getSymbol(), "O")
	TEST_REAL_EQUAL(it->getPosition().x, 3.0)
	TEST_REAL_EQUAL(it->getPosition().y, 4.0)
	TEST_REAL_EQUAL(it->getPosition().z, 5.0)
RESULT

String filename;
NEW_TMP_FILE(filename)
CHECK(XYZFile::write(const System& system))
	System S;
	Molecule* m = new Molecule;
	Atom* a1 = new Atom;
	Atom* a2 = new Atom;
	a1->setElement(PTE[Element::C]);
	a2->setElement(PTE[Element::O]);
	a1->setPosition(Vector3(1.0, 2.0, 3.0));
	a2->setPosition(Vector3(4.0, 5.0, 6.0));
	m->insert(*a1);
	m->insert(*a2);
	S.insert(*m);
	
	XYZFile f(filename, std::ios::out);
	f.write(S);
	f.close();
		
	TEST_FILE(filename.c_str(), "data/XYZFile_test2.xyz")
RESULT


CHECK(XYZFile::XYZFile& operator >> (System& system))
  XYZFile f("data/XYZFile_test1.xyz");
	System S;
	f >> S;
	f.close();

	TEST_EQUAL(S.countAtoms(), 2)
	TEST_EQUAL(S.getName(), "TEST TEST")
	AtomIterator it = S.beginAtom();
	TEST_EQUAL(it->getElement().getSymbol(), "H")
	TEST_REAL_EQUAL(it->getPosition().x, 0.0)
	TEST_REAL_EQUAL(it->getPosition().y, 1.0)
	TEST_REAL_EQUAL(it->getPosition().z, 2.0)
	it++;
	TEST_EQUAL(it->getElement().getSymbol(), "O")
	TEST_REAL_EQUAL(it->getPosition().x, 3.0)
	TEST_REAL_EQUAL(it->getPosition().y, 4.0)
	TEST_REAL_EQUAL(it->getPosition().z, 5.0)
RESULT

CHECK(XYZFile::XYZFile& operator << (const System& system))
	System S;
	Molecule* m = new Molecule;
	Atom* a1 = new Atom;
	Atom* a2 = new Atom;
	a1->setElement(PTE[Element::C]);
	a2->setElement(PTE[Element::O]);
	a1->setPosition(Vector3(1.0, 2.0, 3.0));
	a2->setPosition(Vector3(4.0, 5.0, 6.0));
	m->insert(*a1);
	m->insert(*a2);
	S.insert(*m);
	
	XYZFile f(filename, std::ios::out);
	f << S;
	f.close();
		
	TEST_FILE(filename.c_str(), "data/XYZFile_test2.xyz")
RESULT

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
