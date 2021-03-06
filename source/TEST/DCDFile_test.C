// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// $Id: DCDFile_test.C,v 1.27 2004/11/07 19:54:59 oliver Exp $
//

#include <BALL/CONCEPT/classTest.h>

///////////////////////////
#include <BALL/FORMAT/DCDFile.h>
#include <BALL/FORMAT/PDBFile.h>
#include <BALL/MOLMEC/COMMON/snapShotManager.h>
#include <BALL/MOLMEC/COMMON/snapShot.h>
#include <BALL/MOLMEC/AMBER/amber.h>
///////////////////////////

START_TEST(DCDFile, "$Id: DCDFile_test.C,v 1.27 2004/11/07 19:54:59 oliver Exp $")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace BALL;

String dcd_test_file("data/DCD_test.dcd");

DCDFile* p = 0;
CHECK(DCDFile() throw())
  p = new DCDFile;
	TEST_NOT_EQUAL(p, 0)
RESULT


CHECK(~DCDFile() throw())
	delete p;
RESULT


CHECK(DCDFile(const DCDFile& file) throw(Exception::FileNotFound))
  DCDFile first_file(dcd_test_file);
	DCDFile second_file(first_file);
	bool test = (first_file == second_file);
	TEST_EQUAL(test, true)
RESULT


CHECK(DCDFile(const String& name, File::OpenMode open_mode = std::ios::in) throw())
  DCDFile test_file(dcd_test_file, std::ios::in);
	TEST_EQUAL(test_file.isOpen(), true)
	TEST_EQUAL(test_file.getOpenMode(), std::ios::binary | std::ios::in)
	test_file.close();
	TEST_EQUAL(test_file.isOpen(), false)

	String filename;
	NEW_TMP_FILE(filename)
	DCDFile test_file2(filename, std::ios::out);
	TEST_EQUAL(test_file2.isOpen(), true)
	TEST_EQUAL(test_file2.getOpenMode(), std::ios::binary | std::ios::out)
	test_file2.close();
	TEST_EQUAL(test_file2.isOpen(), false)
RESULT


CHECK(const DCDFile& operator = (const DCDFile& file) throw())
  DCDFile one(dcd_test_file, std::ios::in);
	DCDFile two;
	two = one;
	bool test = (one == two);
	TEST_EQUAL(test, true)
RESULT


CHECK(void clear() throw())
	DCDFile file(dcd_test_file, std::ios::in);
	file.clear();
	DCDFile empty;
	bool test = (file == empty);
	TEST_EQUAL(test, true)
RESULT


CHECK(bool operator == (const DCDFile& file) const throw())
  DCDFile one(dcd_test_file, std::ios::in);
  DCDFile two(dcd_test_file, std::ios::in);
	bool test = (one == two);
	TEST_EQUAL(test, true)
	DCDFile three;
	test = (one == three);
	TEST_EQUAL(test, false)
RESULT

String filename;
System system;
Size nr_of_atoms;

CHECK([EXTRA] full test writing)
	PDBFile pfile("data/DCDFile_test.pdb");
	pfile.read(system);
	nr_of_atoms = system.countAtoms();
	TEST_EQUAL(nr_of_atoms, 892)
	TEST_EQUAL(system.getAtom(0)->getPosition(), Vector3(11.936, 104.294, 10.149))
	AmberFF amberFF;
	NEW_TMP_FILE(filename);
	DCDFile dcd(filename, std::ios::out);
	dcd.enableVelocityStorage();
	TEST_EQUAL(dcd.isAccessible(), true)
	Options options;
	SnapShotManager ssm(&system, &amberFF, options, &dcd);
	ssm.takeSnapShot();
	system.getAtom(0)->setPosition(Vector3(1,2,1111));
	system.getAtom(0)->setForce(Vector3(3,4,5));
	system.getAtom(0)->setVelocity(Vector3(6,7,8));
	ssm.takeSnapShot();
	ssm.flushToDisk();
	dcd.close();
	TEST_EQUAL(dcd.getNumberOfSnapShots(), 2)
	TEST_NOT_EQUAL(dcd.getSize(), 0) 
	system.getAtom(0)->setVelocity(Vector3(1,2,3));
RESULT


SnapShot ss;
CHECK(bool read(SnapShot& snapshot) throw())
	DCDFile dcd(filename);
	TEST_EQUAL(dcd.getNumberOfSnapShots(), 2)
	bool result = dcd.read(ss);
	TEST_EQUAL(result, true)
	ss.applySnapShot(system);
	TEST_EQUAL(ss.getNumberOfAtoms(), nr_of_atoms)
	TEST_EQUAL(ss.getAtomPositions()[0], Vector3(11.936, 104.294, 10.149))
	TEST_EQUAL(system.getAtom(0)->getPosition(), Vector3(11.936, 104.294, 10.149))
	TEST_EQUAL(system.getAtom(0)->getVelocity(), Vector3(0,0,0))
	result = dcd.read(ss);
	TEST_EQUAL(result, true)
	TEST_EQUAL(ss.getNumberOfAtoms(), nr_of_atoms)
	ss.applySnapShot(system);
	TEST_EQUAL(system.getAtom(0)->getPosition(), Vector3(1,2,1111))
	TEST_EQUAL(system.getAtom(0)->getForce(), Vector3(3,4,5))
	TEST_EQUAL(system.getAtom(0)->getVelocity(), Vector3(6,7,8))
	result = dcd.read(ss);
	TEST_EQUAL(result, false)
RESULT


CHECK(bool readHeader() throw())
  DCDFile one(dcd_test_file, std::ios::in);
	bool test = one.readHeader();
	TEST_EQUAL(test, false)
	one.close();
	DCDFile two("data/INIFile_test.ini", std::ios::in);
	test = two.readHeader();
	TEST_EQUAL(test, false)

	DCDFile three(filename);
	test = three.readHeader();
	TEST_EQUAL(test, false)
RESULT


CHECK(bool writeHeader() throw())
	String temporary;
	NEW_TMP_FILE(temporary)
	DCDFile one(temporary, std::ios::out);
	bool test = one.writeHeader();
	TEST_EQUAL(test, true);
	one.close();
	TEST_NOT_EQUAL(one.getSize(), 0)  
	DCDFile two(temporary, std::ios::in);
	test = two.readHeader();
	TEST_EQUAL(test, true);  // got false
RESULT


CHECK(bool append(const SnapShot& snapshot) throw())
	System system;
	Molecule mol;
	Atom atom;

	atom.setPosition(Vector3(1.0, 2.0, 3.0));

	mol.insert(atom);
	system.insert(mol);

	SnapShot snap;
	snap.takeSnapShot(system);

	String temporary;
	NEW_TMP_FILE(temporary)
	DCDFile one(temporary, std::ios::out);
	for (Size i = 0; i < 100; i++) one.append(snap);
	one.close();
	DCDFile two(temporary, std::ios::in);
	SnapShot snap2;
	two.read(snap2);
	// ?????
	// DCDFile lacks something like updateHeader(). As long as such a
	// function doesn't exist, appending cannot work.
	// bool test = (snap == snap2);
	// TEST_EQUAL(test, true);
RESULT


CHECK(bool flushToDisk(const std::vector<SnapShot>& buffer) throw(File::CannotWrite))
	vector<SnapShot> v;
	v.push_back(ss);
	TEST_EQUAL(ss.getNumberOfAtoms(), 892)
	String temporary;
	NEW_TMP_FILE(temporary)
	DCDFile dcd(temporary, std::ios::out);
	TEST_EQUAL(dcd.isOpen(), true)
	TEST_EQUAL(dcd.isWritable(), true)
	TEST_EQUAL(dcd.getOpenMode(), std::ios::binary | std::ios::out)
	TEST_EQUAL(dcd.isAccessible(), true)
	TEST_EQUAL(dcd.getNumberOfSnapShots(), 0)
	bool result = dcd.flushToDisk(v);
	TEST_EQUAL(result, true)
	TEST_EQUAL(dcd.getNumberOfSnapShots(), 1)

	result = dcd.flushToDisk(v);
	TEST_EQUAL(result, true)
	TEST_EQUAL(dcd.getNumberOfSnapShots(), 2)
	dcd.close();
	TEST_NOT_EQUAL(dcd.getSize(), 0)   // got 0

	DCDFile empty;
	TEST_EXCEPTION(File::CannotWrite, empty.flushToDisk(v))
RESULT

CHECK(BALL_CREATE(DCDFile))
	DCDFile dcd(filename);
	DCDFile* dcd_ptr = (DCDFile*)dcd.create();
	TEST_NOT_EQUAL(dcd_ptr, 0)
	TEST_EQUAL(dcd_ptr->getName(), filename)
	delete dcd_ptr;
RESULT

CHECK(bool init() throw())
	// init() is called every time an object is constructed, so we won't have a
	// dedicated test here. There is no reasonable way to test it, anyway.
	DCDFile d;
	d.init();
RESULT

CHECK(bool isSwappingBytes() const throw())
	DCDFile d;
	d.isSwappingBytes();
RESULT

CHECK(bool open(const String& name, File::OpenMode open_mode = std::ios::in) throw(Exception::FileNotFound))
	DCDFile dcd;
	TEST_EXCEPTION(Exception::FileNotFound, dcd.open("aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa", std::ios::in))
RESULT

CHECK(bool seekAndWriteHeader() throw())
	DCDFile d;
	TEST_EQUAL(d.seekAndWriteHeader(), false)
RESULT

CHECK(bool hasVelocities() const throw())
	DCDFile d;
	TEST_EQUAL(d.hasVelocities(), false)
RESULT

CHECK(void disableVelocityStorage() throw())
	DCDFile d;
	d.enableVelocityStorage();
	d.disableVelocityStorage();
	TEST_EQUAL(d.hasVelocities(), false)
RESULT

CHECK(void enableVelocityStorage() throw())
	DCDFile d;
	d.disableVelocityStorage();
	d.enableVelocityStorage();
	TEST_EQUAL(d.hasVelocities(), true)
RESULT

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
