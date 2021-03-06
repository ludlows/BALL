// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// $Id: HashMap_test.C,v 1.16 2003/06/26 10:41:05 amoll Exp $

#include <BALL/CONCEPT/classTest.h>

///////////////////////////
#include <BALL/DATATYPE/hashMap.h>
#include <BALL/CONCEPT/visitor.h>
#include "ItemCollector.h"
///////////////////////////

using namespace BALL;
using namespace std;

class MyVisitor 
	: public Visitor <HashMap<int, int> >
{
	public:
	MyVisitor()
		: value_sum(0),
			key_sum(0)
	{
	}

	int value_sum;
	int key_sum;

	void visit(HashMap<int, int>& v)
	{
		HashMap<int, int>::Iterator it = v.begin();
		for (; it != v.end(); ++it)
		{
			key_sum += (*it).first;
			value_sum += (*it).second;
		}
	}
};

START_TEST(HashMap, "$Id: HashMap_test.C,v 1.16 2003/06/26 10:41:05 amoll Exp $")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

HashMap<int, int>* map_ptr;
CHECK(HashMap(Size initial_capacity = INITIAL_CAPACITY, Size number_of_buckets = INITIAL_NUMBER_OF_BUCKETS) throw())
	map_ptr = new HashMap<int, int>;
	TEST_NOT_EQUAL(map_ptr, 0)
RESULT

CHECK(~HashMap() throw())
	delete map_ptr;
RESULT

CHECK(Size size() const throw())
	HashMap<int, int> hm;
	TEST_EQUAL(hm.size(), 0)
RESULT

CHECK(Size getSize() const throw())
	HashMap<int, int> hm;
	TEST_EQUAL(hm.getSize(), 0)
RESULT

CHECK(Size getBucketSize() const throw())
	HashMap<int, int> hm;
	TEST_EQUAL(hm.getBucketSize(), 3)
RESULT

CHECK(Size getCapacity() const throw())
	HashMap<int, int> hm;
	TEST_EQUAL(hm.getCapacity(), 4)
RESULT

//CHECK(std::pair<Iterator, bool> insert(const ValueType& entry) throw())
CHECK([EXTRA] otherwise_there_is_a_problem_here)
	HashMap<int, int> hm;
	TEST_EQUAL(hm.getSize(), 0)
	hm.insert(HashMap<int, int>::ValueType(0, 0));
	TEST_EQUAL(hm.getSize(), 1)
	hm.insert(HashMap<int, int>::ValueType(0, 1));
	TEST_EQUAL(hm.getSize(), 1)
	hm.insert(HashMap<int, int>::ValueType(1, 0));
	TEST_EQUAL(hm.getSize(), 2)
	hm.insert(HashMap<int, int>::ValueType(1, 1));
	TEST_EQUAL(hm.getSize(), 2)
RESULT

CHECK(Iterator find(const Key& key) throw())
	HashMap<int, int> hm;
	TEST_EQUAL((hm.find(0) == hm.end()), true)
	TEST_EQUAL((hm.find(1) == hm.end()), true)
	TEST_EQUAL((hm.find(2) == hm.end()), true)
	TEST_EQUAL((hm.find(-2) == hm.end()), true)
	hm.insert(HashMap<int, int>::ValueType(0, 0));
	TEST_EQUAL((hm.find(0) == hm.end()), false)
	TEST_EQUAL((hm.find(1) == hm.end()), true)
	TEST_EQUAL((hm.find(2) == hm.end()), true)
	TEST_EQUAL((hm.find(-2) == hm.end()), true)
	hm.insert(HashMap<int, int>::ValueType(1, 1));
	TEST_EQUAL((hm.find(0) == hm.end()), false)
	TEST_EQUAL((hm.find(1) == hm.end()), false)
	TEST_EQUAL((hm.find(2) == hm.end()), true)
	TEST_EQUAL((hm.find(-2) == hm.end()), true)
RESULT

CHECK(ConstIterator find(const Key& key) const throw())
	HashMap<int, int> h_mutable;
	const HashMap<int, int>& hm = const_cast<const HashMap<int, int>&>(h_mutable);
	TEST_EQUAL((hm.find(0) == hm.end()), true)
	TEST_EQUAL((hm.find(1) == hm.end()), true)
	TEST_EQUAL((hm.find(2) == hm.end()), true)
	TEST_EQUAL((hm.find(-2) == hm.end()), true)
	h_mutable.insert(HashMap<int, int>::ValueType(0, 0));
	TEST_EQUAL((hm.find(0) == hm.end()), false)
	TEST_EQUAL((hm.find(1) == hm.end()), true)
	TEST_EQUAL((hm.find(2) == hm.end()), true)
	TEST_EQUAL((hm.find(-2) == hm.end()), true)
	h_mutable.insert(HashMap<int, int>::ValueType(1, 1));
	TEST_EQUAL((hm.find(0) == hm.end()), false)
	TEST_EQUAL((hm.find(1) == hm.end()), false)
	TEST_EQUAL((hm.find(2) == hm.end()), true)
	TEST_EQUAL((hm.find(-2) == hm.end()), true)
RESULT

CHECK(Size erase(const Key& key) throw())
	HashMap<int, int> hm;
	TEST_EQUAL(hm.getSize(), 0)
	hm.insert(HashMap<int, int>::ValueType(0, 0));
	TEST_EQUAL(hm.getSize(), 1)
	hm.insert(HashMap<int, int>::ValueType(1, 1));
	TEST_EQUAL(hm.getSize(), 2)
	hm.insert(HashMap<int, int>::ValueType(2, 2));
	TEST_EQUAL(hm.find(1) == hm.end(), false)
	hm.erase(1);
	TEST_EQUAL(hm.find(1) == hm.end(), true)
	TEST_EQUAL(hm.find(0) == hm.end(), false)
	TEST_EQUAL(hm.find(2) == hm.end(), false)
	TEST_EQUAL(hm.find(5) == hm.end(), true)
	TEST_EQUAL(hm.getSize(), 2)
	hm.erase(5);
	TEST_EQUAL(hm.find(1) == hm.end(), true)
	TEST_EQUAL(hm.find(0) == hm.end(), false)
	TEST_EQUAL(hm.find(2) == hm.end(), false)
	TEST_EQUAL(hm.find(5) == hm.end(), true)
	TEST_EQUAL(hm.getSize(), 2)
RESULT

CHECK(void erase(Iterator first, Iterator last) throw(Exception::IncompatibleIterators))
	HashMap<int, int> hm;
	hm.insert(HashMap<int, int>::ValueType(0, 0));
	hm.insert(HashMap<int, int>::ValueType(1, 1));
	hm.insert(HashMap<int, int>::ValueType(2, 2));
	hm.insert(HashMap<int, int>::ValueType(3, 3));

	HashMap<int, int>::Iterator it1 = hm.begin();
	HashMap<int, int>::Iterator it2 = hm.begin();
	++it2;
	++it2;
	++it1;

	hm.erase(it1, it2);
	TEST_EQUAL(hm.has(0), false)
	TEST_EQUAL(hm.has(2), true)
	TEST_EQUAL(hm.has(3), true)
	TEST_EQUAL(hm.has(1), true)
	TEST_EQUAL(hm.getSize(), 3)

	hm.erase(hm.begin(), hm.end());
	TEST_EQUAL(hm.getSize(), 0)

	hm.insert(HashMap<int, int>::ValueType(0, 0));
	hm.insert(HashMap<int, int>::ValueType(1, 1));
	hm.insert(HashMap<int, int>::ValueType(2, 2));
	hm.insert(HashMap<int, int>::ValueType(3, 3));
	TEST_EQUAL(hm.getSize(), 4)
	it1 = hm.begin();
	++it1;
	++it1;
	hm.erase(it1, hm.end());
	TEST_EQUAL(hm.getSize(), 2)
	TEST_EQUAL(hm.has(0), true)
	TEST_EQUAL(hm.has(1), false)
	TEST_EQUAL(hm.has(2), false)
	TEST_EQUAL(hm.has(3), true)

	HashMap<int, int> hm2;
	TEST_EXCEPTION(Exception::IncompatibleIterators, hm.erase(hm.begin(), hm2.end()))
	TEST_EXCEPTION(Exception::IncompatibleIterators, hm.erase(hm2.begin(), hm.end()))
RESULT

CHECK(void erase(Iterator pos) throw(Exception::IncompatibleIterators, Exception::InvalidIterator))
	HashMap<int, int> hm;
	hm.insert(HashMap<int, int>::ValueType(0, 0));
	hm.insert(HashMap<int, int>::ValueType(1, 1));
	hm.insert(HashMap<int, int>::ValueType(2, 2));
	hm.insert(HashMap<int, int>::ValueType(3, 3));
	HashMap<int, int>::Iterator it1 = hm.begin();
	++it1;
	TEST_EQUAL(it1->first, 0)
	hm.erase(it1);
	TEST_EQUAL(hm.has(0), false)
	TEST_EQUAL(hm.has(1), true)
	TEST_EQUAL(hm.has(2), true)
	TEST_EQUAL(hm.has(3), true)
	TEST_EQUAL(hm.getSize(), 3)
	hm.erase(hm.end());
	TEST_EQUAL(hm.getSize(), 3)
	hm.erase(hm.begin());
	TEST_EQUAL(hm.getSize(), 2)
	hm.erase(hm.begin());
	TEST_EQUAL(hm.getSize(), 1)
	hm.erase(hm.begin());
	TEST_EQUAL(hm.getSize(), 0)
	hm.erase(hm.begin());
	TEST_EQUAL(hm.getSize(), 0)

	HashMap<int, int> hm2;
	TEST_EXCEPTION(Exception::IncompatibleIterators, hm2.erase(hm.begin()))
RESULT

CHECK(T& operator [] (const Key& key) throw())
	HashMap<int, int> hm;
	hm[0] = 0;
	hm[0] = 1;
	hm[1] = 2;
	hm[2] = 4;
	hm[3] = 8;
	hm[4] = 16;
	hm[5] = 32;
	TEST_EQUAL(hm.getSize(), 6)
	TEST_EQUAL(hm[0], 1)
	TEST_EQUAL(hm[1], 2)
	TEST_EQUAL(hm[2], 4)
	TEST_EQUAL(hm[3], 8)
	TEST_EQUAL(hm[4], 16)
	TEST_EQUAL(hm[5], 32)
RESULT

CHECK(const T& operator [] (const Key& key) const throw(typename HashMap<Key, T>::IllegalKey))
	HashMap<int, int> hm;
	hm[0] = 0;
	hm[0] = 1;
	hm[1] = 2;
	hm[2] = 4;
	hm[3] = 8;
	hm[4] = 16;
	hm[5] = 32;
	const HashMap<int, int>& const_map = const_cast<const HashMap<int, int>&>(hm);
	TEST_EQUAL(const_map.getSize(), 6)
	TEST_EQUAL(const_map[0], 1)
	TEST_EQUAL(const_map[1], 2)
	TEST_EQUAL(const_map[2], 4)
	TEST_EQUAL(const_map[3], 8)
	TEST_EQUAL(const_map[4], 16)
	TEST_EQUAL(const_map[5], 32)
	typedef HashMap<int,int> MyHashMap; // otherwise next line wont work
	TEST_EXCEPTION(MyHashMap::IllegalKey, const_map[6])
RESULT

CHECK(HashMap(const HashMap& hash_map) throw())
	HashMap<int, int> hm;
	hm.insert(HashMap<int, int>::ValueType(0, 0));
	hm.insert(HashMap<int, int>::ValueType(1, 1));
	hm.insert(HashMap<int, int>::ValueType(2, 2));
	hm.insert(HashMap<int, int>::ValueType(3, 3));

	HashMap<int, int> hm2(hm);
	TEST_EQUAL(hm2[0], 0)
	TEST_EQUAL(hm2[1], 1)
	TEST_EQUAL(hm2[2], 2)
	TEST_EQUAL(hm2[3], 3)
RESULT

CHECK(void clear() throw())
	HashMap<int, int> hm;
	hm.insert(HashMap<int, int>::ValueType(0, 0));
	hm.clear();
	TEST_EQUAL(hm.getSize(), 0)
	TEST_EQUAL(hm.getCapacity(), 4)
	TEST_EQUAL(hm.getBucketSize(), 3)
RESULT

CHECK(void destroy() throw())
	HashMap<int, int> hm;
	hm.insert(HashMap<int, int>::ValueType(0, 0));
	hm.destroy();
	TEST_EQUAL(hm.getSize(), 0)
	TEST_EQUAL(hm.getCapacity(), 4)
	TEST_EQUAL(hm.getBucketSize(), 3)
RESULT

CHECK(void set(const HashMap& hash_map) throw())
	HashMap<int, int> hm1;
	hm1.insert(pair<int, int>(12, 34));
	hm1.insert(pair<int, int>(44, 55));
	HashMap<int, int> hm2;
	TEST_EQUAL(hm1.getSize(), 2)
	TEST_EQUAL(hm2.getSize(), 0)
	TEST_EQUAL(hm1.has(12), true)
	TEST_EQUAL(hm1.has(44), true)
	hm2.set(hm1);
	TEST_EQUAL(hm2.getSize(), 2)
	TEST_EQUAL(hm2.has(12), true)
	TEST_EQUAL(hm2.has(44), true)
	TEST_EQUAL(hm1.getSize(), 2)
	TEST_EQUAL(hm1.has(12), true)
	TEST_EQUAL(hm1.has(44), true)
RESULT

CHECK(void get(HashMap& hash_map) const throw())
	HashMap<int, int> hm, hm2;
	hm.insert(HashMap<int, int>::ValueType(0, 0));
	hm.insert(HashMap<int, int>::ValueType(1, 1));

	hm.get(hm2);
	TEST_EQUAL(hm2[0], 0)
	TEST_EQUAL(hm2[1], 1)
RESULT

CHECK(void swap(HashMap& hash_map) throw())
	HashMap<int, int> hm, hm2;
	hm.insert(HashMap<int, int>::ValueType(0, 0));
	hm.insert(HashMap<int, int>::ValueType(1, 1));

	hm2.insert(HashMap<int, int>::ValueType(10, 10));
	hm2.insert(HashMap<int, int>::ValueType(11, 11));

	TEST_EQUAL(hm.getSize(), 2)
	TEST_EQUAL(hm2.getSize(), 2)
	TEST_EQUAL(hm2[10], 10)
	TEST_EQUAL(hm2[11], 11)
	TEST_EQUAL(hm[0], 0)
	TEST_EQUAL(hm[1], 1)
	hm.swap(hm2);
	TEST_EQUAL(hm.getSize(), 2)
	TEST_EQUAL(hm2.getSize(), 2)
	TEST_EQUAL(hm[10], 10)
	TEST_EQUAL(hm[11], 11)
	TEST_EQUAL(hm2[0], 0)
	TEST_EQUAL(hm2[1], 1)
RESULT

CHECK(const HashMap& operator = (const HashMap& hash_map) throw())
	HashMap<int, int> hm, hm2;
	hm.insert(HashMap<int, int>::ValueType(0, 0));
	hm.insert(HashMap<int, int>::ValueType(1, 1));

	hm2 = hm;
	TEST_EQUAL(hm2[0], 0)
	TEST_EQUAL(hm2[1], 1)
RESULT

CHECK(void host(Visitor<HashMap<Key, T> >& visitor) throw())
	HashMap<int, int> hm;
	hm.insert(HashMap<int, int>::ValueType(100, 200));
	hm.insert(HashMap<int, int>::ValueType(11, 22));

	MyVisitor mv;
	hm.host(mv);
	TEST_EQUAL(mv.key_sum, 111)
	TEST_EQUAL(mv.value_sum, 222)
RESULT

CHECK(bool has(const Key& key) const throw())
	HashMap<int, int> hm;
	hm.insert(HashMap<int, int>::ValueType(0, 0));
	hm.insert(HashMap<int, int>::ValueType(1, 1));
	TEST_EQUAL(hm.has(0), true)
	TEST_EQUAL(hm.has(1), true)
	TEST_EQUAL(hm.has(2), false)
RESULT

CHECK(bool isEmpty() const throw())
	HashMap<int, int> hm;
	TEST_EQUAL(hm.isEmpty(), true)
	hm.insert(HashMap<int, int>::ValueType(0, 0));
	TEST_EQUAL(hm.isEmpty(), false)
	hm.clear();
	TEST_EQUAL(hm.isEmpty(), true)
RESULT

CHECK(bool operator == (const HashMap& hash_map) const throw())
	HashMap<int, int> hm;
	HashMap<int, int> hm2;
	TEST_EQUAL(hm == hm2, true)

	hm.insert(HashMap<int, int>::ValueType(0, 0));
	hm.insert(HashMap<int, int>::ValueType(1, 1));

	hm2.insert(HashMap<int, int>::ValueType(0, 0));
	hm2.insert(HashMap<int, int>::ValueType(1, 1));
	TEST_EQUAL(hm == hm2, true)

	hm2.insert(HashMap<int, int>::ValueType(1, 2));
	TEST_EQUAL(hm[1], 1)
	TEST_EQUAL(hm2[1], 2)
	TEST_EQUAL(hm == hm2, false)

	hm2.insert(HashMap<int, int>::ValueType(1, 1));
	hm2.insert(HashMap<int, int>::ValueType(2, 1));
	TEST_EQUAL(hm == hm2, false)
RESULT


CHECK(bool operator != (const HashMap& hash_map) const throw())
	HashMap<int, int> hm;
	HashMap<int, int> hm2;

	TEST_EQUAL(hm != hm2, false)
	hm.insert(HashMap<int, int>::ValueType(0, 0));
	hm.insert(HashMap<int, int>::ValueType(1, 1));

	hm2.insert(HashMap<int, int>::ValueType(0, 0));
	hm2.insert(HashMap<int, int>::ValueType(1, 1));

	TEST_EQUAL(hm != hm2, false)
	hm2.insert(HashMap<int, int>::ValueType(2, 1));
	TEST_EQUAL(hm != hm2, true)
RESULT


CHECK(void dump(::std::ostream& s = ::std::cout, Size depth = 0) const throw())
	HashMap<int, int> hm;
	hm.insert(HashMap<int, int>::ValueType(0, 0));
	hm.insert(HashMap<int, int>::ValueType(1, 1));

  String filename;
	NEW_TMP_FILE(filename)
	std::ofstream outfile(filename.c_str(), std::ios::out);
	hm.dump(outfile);
	outfile.close();
	TEST_FILE_REGEXP(filename.c_str(), "data/HashMap_test.txt")
RESULT


CHECK(bool isValid() const throw())
	HashMap<int, int> hm;
	TEST_EQUAL(hm.isValid(), true)
RESULT


CHECK(bool apply(UnaryProcessor<ValueType>& processor) throw())
	HashMap<int, int> hm;
	hm.insert(HashMap<int, int>::ValueType(0, 0));
	hm.insert(HashMap<int, int>::ValueType(1, 1));
	ItemCollector<pair<int,int> > myproc;
	TEST_EQUAL(hm.apply(myproc), true)
	TEST_EQUAL(myproc.getSize(), 2)
	TEST_EQUAL(myproc.getPointer()->first, 0) myproc.forward();
	TEST_EQUAL(myproc.getPointer()->first, 1) myproc.forward();
RESULT


CHECK(BALL_CREATE(HashMap))
  // ???
RESULT

CHECK(ConstIterator begin() const throw())
	HashMap<int, int> hm;
	hm.insert(HashMap<int, int>::ValueType(0, 0));
	HashMap<int, int>::ConstIterator it = hm.begin();
	TEST_EQUAL(it->first == 0, true)
	it++;
	TEST_EQUAL(it == hm.end(), true)
RESULT

CHECK(ConstIterator end() const throw())
	HashMap<int, int> hm;
	TEST_EQUAL(hm.begin() ==  hm.end(), true)
RESULT

CHECK(HashMap* getContainer() throw())
  // ???
RESULT

CHECK(IllegalKey(const char* file, int line))
	HashMap<int, int>::IllegalKey ik(__FILE__, __LINE__);
RESULT

CHECK(Iterator end() throw())
	HashMap<int, int> hm;
	HashMap<int, int>::Iterator it1 = hm.begin();
	HashMap<int, int>::Iterator it2 = hm.end();
	TEST_EQUAL(it1 == it2, true)
RESULT

CHECK(Iterator insert(Iterator pos, const ValueType& entry) throw())
  // ???
RESULT

HashMap<int, int> hm;
hm.insert(HashMap<int, int>::ValueType(0, 0));
hm.insert(HashMap<int, int>::ValueType(1, 1));
HashMap<int, int>::Iterator it = hm.begin();


// IteratorTraits_ tests
CHECK(IteratorPosition& getPosition() throw())
  // ???
RESULT

CHECK(IteratorTraits_() throw())
  // ???
RESULT

CHECK(IteratorTraits_(const HashMap& hash_map) throw())
  // ???
RESULT

CHECK(IteratorTraits_(const IteratorTraits_& traits) throw())
  // ???
RESULT

CHECK(ValueType& getData() throw())
  // ???
RESULT

CHECK(bool isBegin() const throw())
  // ???
RESULT

CHECK(bool isEnd() const throw())
  // ???
RESULT

CHECK(bool isSingular() const throw())
  // ???
RESULT

CHECK(bool operator != (const IteratorTraits_& traits) const throw())
  // ???
RESULT

CHECK(bool operator == (const IteratorTraits_& traits) const throw())
  // ???
RESULT

CHECK(const HashMap* getContainer() const throw())
  // ???
RESULT

CHECK(const IteratorPosition& getPosition() const throw())
  // ???
RESULT

CHECK(const IteratorTraits_& operator = (const IteratorTraits_& traits) throw())
  // ???
RESULT

CHECK(const ValueType& getData() const throw())
  // ???
RESULT

CHECK(friend Iterator begin() throw())
  // ???
RESULT

CHECK((std::pair<Iterator, bool> insert(const ValueType& entry) throw()))
  // ???
RESULT

CHECK(void forward() throw())
  // ???
RESULT

CHECK(void invalidate() throw())
  // ???
RESULT

CHECK(void toBegin() throw())
  // ???
RESULT

CHECK(void toEnd() throw())
  // ???
RESULT

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
