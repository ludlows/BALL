// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// $Id: LogStream_test.C,v 1.22.4.4 2005/10/21 13:01:01 amoll Exp $
//

#include <BALL/CONCEPT/classTest.h>

///////////////////////////

#include <BALL/COMMON/logStream.h>
#ifdef BALL_HAS_SYS_TIME_H
#	include <sys/time.h>
#endif
#include <BALL/MATHS/common.h>

///////////////////////////

START_TEST(LogStream, "$Id: LogStream_test.C,v 1.22.4.4 2005/10/21 13:01:01 amoll Exp $")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////


using namespace BALL;
using namespace std;

class TestTarget
	:	public LogStreamNotifier
{
	public:
	virtual void logNotify()
	{
		notified = true;
		return;
	}
	bool notified;
};

String filename;

LogStream* l1 = 0;

CHECK(LogStream(LogStreamBuf*, bool, bool))
	l1 = new LogStream(new LogStreamBuf);
	TEST_NOT_EQUAL(l1, 0)
RESULT

CHECK(~LogStream())
	delete l1;
RESULT

CHECK(LogStream(LogStreamBuf* buf))
	LogStream* l1;
	l1 = new LogStream((LogStreamBuf*)0);
	TEST_NOT_EQUAL(l1, 0)
	delete l1;

	LogStream* l2;
	LogStreamBuf* lb2 = new LogStreamBuf;
	l2 = new LogStream(lb2);
	TEST_NOT_EQUAL(l2, 0)
	delete l2;
RESULT

CHECK(rdbuf())
	LogStream l1(new LogStreamBuf);
	TEST_NOT_EQUAL(l1.rdbuf(), 0)
RESULT

CHECK(operator -> ())
	LogStream l1(new LogStreamBuf);
	l1->sync();
RESULT

CHECK(setLevel(int level))
	LogStream l1(new LogStreamBuf);
	l1 << "TEST1" << endl;
	l1.setLevel(99);
	l1 << "TEST2" << endl;
	TEST_EQUAL(l1.getLineLevel(0), 0)
	TEST_EQUAL(l1.getLineLevel(1), 99)
RESULT

CHECK(l1.getLevel())
	LogStream l1(new LogStreamBuf);
	TEST_EQUAL(l1.getLevel(), 0)
	l1.setLevel(99);
	TEST_EQUAL(l1.getLevel(), 99)
RESULT

CHECK(level(int n))
	LogStream l1(new LogStreamBuf);
	l1.level(99) << "TEST" <<endl;
	TEST_EQUAL(l1.getNumberOfLines(), 1)
	TEST_EQUAL(l1.getLineText(0), "TEST")
	TEST_EQUAL(l1.getLineLevel(0), 99)
RESULT

CHECK(info(int n = 0))
	LogStream l1(new LogStreamBuf);
	l1.info() << "TEST" <<endl;
	TEST_EQUAL(l1.getNumberOfLines(), 1)
	TEST_EQUAL(l1.getLineText(0), "TEST")
	TEST_EQUAL(l1.getLineLevel(0), LogStream::INFORMATION_LEVEL)
	l1.info(1) << "TEST2" <<endl;
	TEST_EQUAL(l1.getLineLevel(1), LogStream::INFORMATION_LEVEL + 1)
RESULT

CHECK(error(int n = 0))
	LogStream l1(new LogStreamBuf);
	l1.error() << "TEST" <<endl;
	TEST_EQUAL(l1.getNumberOfLines(), 1)
	TEST_EQUAL(l1.getLineText(0), "TEST")
	TEST_EQUAL(l1.getLineLevel(0), LogStream::ERROR_LEVEL)
	l1.error(1) << "TEST2" <<endl;
	TEST_EQUAL(l1.getLineLevel(1), LogStream::ERROR_LEVEL + 1)
RESULT

CHECK(warn(int n = 0))
	LogStream l1(new LogStreamBuf);
	l1.warn() << "TEST" <<endl;
	TEST_EQUAL(l1.getNumberOfLines(), 1)
	TEST_EQUAL(l1.getLineText(0), "TEST")
	TEST_EQUAL(l1.getLineLevel(0), LogStream::WARNING_LEVEL)
	l1.warn(1) << "TEST2" <<endl;
	TEST_EQUAL(l1.getLineLevel(1), LogStream::WARNING_LEVEL + 1)
RESULT

CHECK(insert(std::ostream& s, int min_level = LogStreamBuf::MIN_LEVEL, int max_level = LogStreamBuf::MAX_LEVEL))
	NEW_TMP_FILE(filename)
	LogStream l1(new LogStreamBuf);
	ofstream s(filename.c_str(), std::ios::out);
	l1.insert(s, 99, 99);
	l1.level(98) << "X" << endl;
	l1.level(99) << "1" << endl;
	l1.info(99)  << "2" << endl;
	l1.level(100)<< "X" << endl;
	TEST_EQUAL(l1.getNumberOfLines(), 4)
	TEST_FILE(filename.c_str(), "data/LogStream_test_general.txt")
RESULT

CHECK(remove(std::ostream& s))
	LogStream l1(new LogStreamBuf);
	ofstream s;
	l1.insert(s);
	l1.remove(s);
	// make sure we can remove it twice without harm
	l1.remove(s);
RESULT

CHECK(void insertNotification())
	LogStream l1(new LogStreamBuf);
	TestTarget target;
	ofstream os;
	target.registerAt(l1);
	target.notified = false;
	TEST_EQUAL(target.notified, false)
	l1 << "test" << std::endl;
	TEST_EQUAL(target.notified, true)
RESULT

CHECK(removeNotification)
	LogStream l1(new LogStreamBuf);
	TestTarget target;
	ofstream os;
	target.registerAt(l1);
	target.unregister();
	target.notified = false;
	TEST_EQUAL(target.notified, false)
	l1 << "test" << endl;
	TEST_EQUAL(target.notified, false)
	// make sure we can remove it twice
	target.unregister();
	l1 << "test" << endl;
	TEST_EQUAL(target.notified, false)
RESULT

CHECK(setMinLevel(const std::ostream& s, int min_level))
	NEW_TMP_FILE(filename)
	LogStream l1(new LogStreamBuf);
	ofstream s(filename.c_str(), std::ios::out);
	l1.insert(s, 0);
	l1.setMinLevel(s, 98);
	l1.info(97) << "X" << endl;
	l1.info(98) << "1" << endl;
	l1.info(99) << "2" << endl;
	TEST_EQUAL(l1.getNumberOfLines(), 3)
	TEST_FILE(filename.c_str(), "data/LogStream_test_general.txt")
RESULT

CHECK(setMaxLevel(const std::ostream& s, int max_level))
	NEW_TMP_FILE(filename)
	LogStream l1(new LogStreamBuf);
	ofstream s(filename.c_str(), std::ios::out);
	l1.insert(s, 0);
	l1.setMaxLevel(s, 98);
	l1.info(97) << "1" << endl;
	l1.info(98) << "2" << endl;
	l1.info(99) << "X" << endl;
	TEST_EQUAL(l1.getNumberOfLines(), 3)
	TEST_FILE(filename.c_str(), "data/LogStream_test_general.txt")
RESULT

CHECK(setPrefix(const std::ostream& s, const string& prefix))
	NEW_TMP_FILE(filename)
	LogStream l1(new LogStreamBuf);
	ofstream s(filename.c_str(), std::ios::out);;
	l1.insert(s, 0);
	l1.setPrefix(s, "%l"); //loglevel
	l1.info(1) << "  1." << endl;
	l1.setPrefix(s, "%y"); //message type ("Error", "Warning", "Information", "-")
	l1.info(2) << "  2." << endl;
	l1.setPrefix(s, "%T"); //time (HH:MM:SS)
	l1.info(3) << "  3." << endl;
	l1.setPrefix(s, "%t"); //time in short format (HH:MM)
	l1.info(4) << "  4." << endl;
	l1.setPrefix(s, "%D"); //date (DD.MM.YYYY)
	l1.info(5) << "  5." << endl;
	l1.setPrefix(s, "%d"); // date in short format (DD.MM.)
	l1.info(6) << "  6." << endl;
	l1.setPrefix(s, "%S"); //time and date (DD.MM.YYYY, HH:MM:SS)
	l1.info(7) << "  7." << endl;
	l1.setPrefix(s, "%s"); //time and date in short format (DD.MM., HH:MM)
	l1.info(8) << "  8." << endl;
	l1.setPrefix(s, "%%"); //percent sign (escape sequence)
	l1.info(9) << "  9." << endl;
	l1.setPrefix(s, ""); //no prefix
	l1.info(10) << " 10." << endl;
	TEST_EQUAL(l1.getNumberOfLines(), 10)
	TEST_FILE_REGEXP(filename.c_str(), "data/LogStream_test_setPrefix.txt")
RESULT

CHECK(clear())
	LogStream l1(new LogStreamBuf);
	l1.error() << "TEST" <<endl;
	l1.clear();
	TEST_EQUAL(l1.getNumberOfLines(), 0)
RESULT

CHECK(getNumberOfLines(int min_level = LogStreamBuf::MIN_LEVEL, int max_level = LogStreamBuf::MAX_LEVEL))
	LogStream l1(new LogStreamBuf);
	TEST_EQUAL(l1.getNumberOfLines(), 0)
	l1.error() << "TEST" <<endl;
	TEST_EQUAL(l1.getNumberOfLines(), 1)
RESULT

CHECK(getLineText(Index index))
	LogStream l1(new LogStreamBuf);
	l1.error() << "TEST" <<endl;
	TEST_EQUAL(l1.getLineText(0), "TEST")
RESULT

CHECK(getLineTime(Index index))
	LogStream l1(new LogStreamBuf);
	Time timer;
  timer = time(NULL);
	l1.error() << "TEST" <<endl;
	TEST_EQUAL(timer, l1.getLineTime(0))
	TEST_EQUAL(Maths::isNear(timer, l1.getLineTime(0), (Time)1), true)
RESULT

CHECK(getLineLevel(Index index))
	LogStream l1(new LogStreamBuf);
	l1.level(99) << "TEST" <<endl;
	TEST_EQUAL(l1.getLineLevel(0), 99)
RESULT

CHECK(filterLines(const int min_level = LogStreamBuf::MIN_LEVEL, const int max_level = LogStreamBuf::MAX_LEVEL,
									const Time earliest = 0, const Time latest = LogStreamBuf::MAX_TIME, const string& s = ""))
	LogStream l1(new LogStreamBuf);
  using std::list;													
	list<int>	liste;

	l1.level(0) << "TEST1" << endl;
	l1.level(1) << "TEST2" << endl;
	l1.level(2) << "XXXXX" << endl;
	Time timer = time(NULL);
	Time x = timer;
	while (x == timer)
	{
		x = time(NULL);
	}
	timer = time(NULL);
	l1.level(3) << "XXXXX" << endl;
	l1.level(4) << "TEST4" << endl;
	TEST_EQUAL(l1.getNumberOfLines(), 5)

	liste = l1.filterLines(1, 2, 0, LogStreamBuf::MAX_TIME, "" );
	TEST_EQUAL(liste.size(), 2)
	TEST_EQUAL(liste.front(), 1)
	liste.pop_front();
	TEST_EQUAL(liste.front(), 2)
	liste.clear();

	liste = l1.filterLines(1, 2, 0, LogStreamBuf::MAX_TIME, "XXXXX" );
	TEST_EQUAL(liste.size(), 1)
	TEST_EQUAL(liste.front(), 2)
	liste.clear();

	liste = l1.filterLines(3, 3, timer, LogStreamBuf::MAX_TIME, "XXXXX");
	TEST_EQUAL(liste.size(), 1)
	TEST_EQUAL(liste.front(), 3)
	liste.clear();	
RESULT

// test for a minimum string length for output
CHECK(Output length)
	LogStream l1(new LogStreamBuf);
	l1	<< "aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa" << endl;
RESULT

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
