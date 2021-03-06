// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// $Id: KernelPredicate_test.C,v 1.4 2002/02/27 12:24:36 sturm Exp $
#include <BALL/CONCEPT/classTest.h>

///////////////////////////
#include <BALL/KERNEL/predicate.h>
#include <BALL/KERNEL/system.h>
#include <BALL/KERNEL/molecule.h>
#include <BALL/KERNEL/atom.h>
#include <BALL/KERNEL/bond.h>
#include <BALL/KERNEL/atomContainer.h>
///////////////////////////

START_TEST(KernelPredicate<T>, "$Id: KernelPredicate_test.C,v 1.4 2002/02/27 12:24:36 sturm Exp $")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace BALL;
	
CHECK(KernelPredicate::operator () (const Composite&) const)
System S;
Molecule M;
Atom A;
M.insert(A);
S.insert(M);

KernelPredicate<System>				system_predicate;
KernelPredicate<Molecule>			molecule_predicate;
KernelPredicate<Atom>					atom_predicate;
KernelPredicate<Bond>					bond_predicate;
KernelPredicate<AtomContainer>	atom_container_predicate;
TEST_EQUAL(system_predicate(S), true)
TEST_EQUAL(system_predicate(M), false)
TEST_EQUAL(system_predicate(A), false)
TEST_EQUAL(molecule_predicate(S), false)
TEST_EQUAL(molecule_predicate(M), true)
TEST_EQUAL(molecule_predicate(A), false)
TEST_EQUAL(atom_predicate(S), false)
TEST_EQUAL(atom_predicate(M), false)
TEST_EQUAL(atom_predicate(A), true)
TEST_EQUAL(bond_predicate(S), false)
TEST_EQUAL(bond_predicate(M), false)
TEST_EQUAL(bond_predicate(A), false)
TEST_EQUAL(atom_container_predicate(S), true)
TEST_EQUAL(atom_container_predicate(M), true)
TEST_EQUAL(atom_container_predicate(A), false)
RESULT

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
