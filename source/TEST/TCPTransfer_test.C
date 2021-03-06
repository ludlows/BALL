// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// $Id: TCPTransfer_test.C,v 1.21.4.4 2005/06/20 19:50:10 oliver Exp $
//

#include <BALL/CONCEPT/classTest.h>

///////////////////////////

#include <fstream>

#include <BALL/SYSTEM/TCPTransfer.h>
#include <BALL/SYSTEM/file.h>

///////////////////////////

using namespace BALL;
using namespace std;

// include a helper class testing for network availability
// (or more precisely: the ability to perform HTTP or FTP transfers)

#include "networkTest.h"

START_TEST(TCPTransfer, "$Id: TCPTransfer_test.C,v 1.21.4.4 2005/06/20 19:50:10 oliver Exp $")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
String filename;

TCPTransfer* tcp_ptr = 0;
CHECK(cstr)
	tcp_ptr = new TCPTransfer;
	TEST_NOT_EQUAL(tcp_ptr, 0)
	ABORT_IF(tcp_ptr == 0)
	TEST_EQUAL(tcp_ptr->getStatusCode(), TCPTransfer::UNINITIALIZED__ERROR)
RESULT

CHECK(~TCPTransfer_test)
	delete tcp_ptr;
RESULT

CHECK(http/no login)
	ABORT_IF(!NetworkTest::test("www.ball-project.org", NetworkTest::HTTP))
	NEW_TMP_FILE(filename)
	std::ofstream os(filename.c_str(), std::ios::out);
	
	TCPTransfer tcp_t(os ,"http://www.bioinf.uni-sb.de/OK/BALL/Downloads/http_test.txt");
	TEST_EQUAL(tcp_t.getHostAddress(), "www.bioinf.uni-sb.de")
	TEST_EQUAL(tcp_t.getFileAddress(), "/OK/BALL/Downloads/http_test.txt")
	TEST_EQUAL(tcp_t.getPort(), 80)
	TEST_EQUAL(tcp_t.getStatusCode(), TCPTransfer::OK)
	TEST_EQUAL(tcp_t.getReceivedBytes(), 3048)
	TEST_EQUAL(tcp_t.getLogin(), "")
	TEST_EQUAL(tcp_t.getPassword(), "")
	TEST_EQUAL(tcp_t.getStream(), &os)
	os.close();
	
	TEST_FILE(filename.c_str(), "data/http_test.txt")
RESULT

CHECK(set(ofstream& file, const String& address))
	ABORT_IF(!NetworkTest::test("www.ball-project.org", NetworkTest::HTTP))
	String filename;
	NEW_TMP_FILE(filename)
	std::ofstream os(filename.c_str(), std::ios::out);
	
	TCPTransfer tcp_t;
	tcp_t.set(os ,"http://www.bioinf.uni-sb.de/OK/BALL/Downloads/http_test.txt");
	TEST_EQUAL(tcp_t.getHostAddress(), "www.bioinf.uni-sb.de")
	TEST_EQUAL(tcp_t.getFileAddress(), "/OK/BALL/Downloads/http_test.txt")
	TEST_EQUAL(tcp_t.getPort(), 80)
	TEST_EQUAL(tcp_t.getStatusCode(), TCPTransfer::OK)
	TEST_EQUAL(tcp_t.getReceivedBytes(), 0)
	TEST_EQUAL(tcp_t.getLogin(), "")
	TEST_EQUAL(tcp_t.getPassword(), "")
	TEST_EQUAL(tcp_t.getStream(), &os)	
	tcp_t.transfer();
	TEST_EQUAL(tcp_t.getReceivedBytes(), 3048)
	os.close();

	TEST_FILE(filename.c_str(), "data/http_test.txt")
RESULT

CHECK(http/login)
	ABORT_IF(!NetworkTest::test("www.zbi.uni-saarland.de", NetworkTest::HTTP))
	NEW_TMP_FILE(filename)
	std::ofstream os(filename.c_str(), std::ios::out);
	
	TCPTransfer tcp_t(os ,"http://BALL-TEST:test@www.zbi.uni-saarland.de/zbi/BALL/test/protected/TCPTransferTest.txt");
	TEST_EQUAL(tcp_t.getHostAddress(), "www.zbi.uni-saarland.de")
	TEST_EQUAL(tcp_t.getFileAddress(), "/zbi/BALL/test/protected/TCPTransferTest.txt")
	TEST_EQUAL(tcp_t.getPort(), 80)
	TEST_EQUAL(tcp_t.getStatusCode(), TCPTransfer::OK)
	TEST_EQUAL(tcp_t.getReceivedBytes(), 3048)
	TEST_EQUAL(tcp_t.getLogin(), "BALL-TEST")
	TEST_EQUAL(tcp_t.getPassword(), "test")	
	os.close();

	TEST_FILE(filename.c_str(), "data/http_test.txt")
RESULT

CHECK(ftp)
	ABORT_IF(!NetworkTest::test("ftp.mpi-sb.mpg.de", NetworkTest::FTP))
	NEW_TMP_FILE(filename);
	std::ofstream os(filename.c_str(), std::ios::out);
	
	TCPTransfer tcp_t(os, "ftp://ftp.mpi-sb.mpg.de/pub/outgoing/BALL/ftp_test.txt");
	TEST_EQUAL(tcp_t.getHostAddress(), "ftp.mpi-sb.mpg.de")
	TEST_EQUAL(tcp_t.getFileAddress(), "/pub/outgoing/BALL/ftp_test.txt")
	TEST_EQUAL(tcp_t.getPort(), 21)
	TEST_EQUAL(tcp_t.getStatusCode(), TCPTransfer::OK)
	TEST_EQUAL(tcp_t.getReceivedBytes(), 2312)
	TEST_EQUAL(tcp_t.getLogin(), "")
	TEST_EQUAL(tcp_t.getPassword(), "")
	TEST_EQUAL(tcp_t.getStream(), &os)

	os.close();

	TEST_FILE(filename.c_str(), "data/ftp_test.txt")
RESULT

CHECK(http/exception)
	NEW_TMP_FILE(filename)
	std::ofstream os(filename.c_str(), std::ios::out);
	TEST_EXCEPTION(TCPTransfer::TransferFailed, TCPTransfer tcp_t(os, "ftp://xcajsjddnnakadnndakndna.de/ffaadad.caadd"))
	os.close();
RESULT

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
