// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// $Id: Directory_test.C,v 1.16.4.5 2005/12/06 12:43:09 amoll Exp $
//

#include <BALL/CONCEPT/classTest.h>

///////////////////////////
#include <BALL/SYSTEM/directory.h>
///////////////////////////

#ifdef BALL_PLATFORM_WINDOWS
#	define chdir(a) _chdir(a)
#	define getcwd(a, b) _getcwd(a, b)
#endif


#ifdef BALL_COMPILER_MSVC
# include<windows.h>
# define sleep(a) Sleep(1000 * a)
#endif

using namespace BALL;

String PS = FileSystem::PATH_SEPARATOR;
String test_dir;

bool cleanup()
{
	if (test_dir == "")
	{
		char buffer[10000];
		::getcwd(buffer, Directory::MAX_PATH_LENGTH);
		test_dir = buffer;
		#ifdef BALL_PLATFORM_WINDOWS
			test_dir += "\\data\\Directory_test";
		#else
			test_dir += "/data/Directory_test";
		#endif

		::chdir(test_dir.c_str());
	}

	if (::chdir(test_dir.c_str()) != 0)
	{
		Log.error() << "Could not chdir to test dir!" << std::endl;
	}

	Directory d(test_dir + PS + "dir_a" + PS + "dir_c");

	d.remove("test1" + PS + "test1");
	d.remove("test1" + PS + "test2");
	d.remove("test1" + PS + "test3");
	d.remove("test1" + PS + "test4");
	d.remove("test1");

	d.remove("test2" + PS + "test1");
	d.remove("test2" + PS + "test2");
	d.remove("test2" + PS + "test3");
	d.remove("test2" + PS + "test4");
	d.remove("test2");
	
	//d.remove("test1" + PS + "test2");
	d.remove("test3" + PS + "test1");
	d.remove("test3" + PS + "test2");
	d.remove("test3" + PS + "test3");
	d.remove("test3" + PS + "test4");
	d.remove("test3");

	if (d.has("test1") ||
			d.has("test2") ||
			d.has("test3"))
	{
				return false;

	}

	return true;
}


START_TEST(Directory, "$Id: Directory_test.C,v 1.16.4.5 2005/12/06 12:43:09 amoll Exp $")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace BALL;

CHECK([EXTRA]prerequisites)
	bool ok = cleanup();
	TEST_EQUAL(ok, true)
RESULT

Directory* dd = 0;
char buffer[10000];

CHECK(Directory())
	dd = new Directory();
	TEST_NOT_EQUAL(dd, 0)
	::getcwd(buffer, Directory::MAX_PATH_LENGTH);
	TEST_EQUAL(dd->getPath(), String(buffer))
	TEST_EQUAL(dd->isCurrent(), true)
RESULT

CHECK(~Directory())
	delete dd;
RESULT

CHECK(bool setCurrent(String directory_path))
	Directory d;
	bool result = d.setCurrent(test_dir);
	TEST_EQUAL(result, true);
	::getcwd(buffer, Directory::MAX_PATH_LENGTH);
	String path  = buffer;
	TEST_EQUAL(path, test_dir)
	result = d.setCurrent("dir_a");
	TEST_EQUAL(result, true);
	::getcwd(buffer, Directory::MAX_PATH_LENGTH);
	path = buffer;
	
	const String PS = FileSystem::PATH_SEPARATOR;
	TEST_EQUAL(path, test_dir + PS+"dir_a")
	result = d.setCurrent("c");
	TEST_EQUAL(result, false);
	::getcwd(buffer, Directory::MAX_PATH_LENGTH);
	path = buffer;
	TEST_EQUAL(path, test_dir + PS+ "dir_a")
RESULT


CHECK(Directory(const String& directory_path, bool set_current = false))
	::chdir(test_dir.c_str());
	Directory d("dir_a");
	TEST_EQUAL(d.isValid(), true)

	::getcwd(buffer, Directory::MAX_PATH_LENGTH);
	String s = String(buffer) + PS + "dir_a";
	TEST_EQUAL(d.getPath(), s)

	Directory d1("dir_a"+PS, true);
	TEST_EQUAL(d1.isCurrent(), true)
	TEST_EQUAL(d1.isValid(), true)
	TEST_EQUAL(d.setCurrent(test_dir), true);

	Directory d2("dir_a" + PS + "dir_c"+PS);
	TEST_EQUAL(d2.isValid(), true)

	Directory d3(test_dir);
	TEST_EQUAL(d3.isValid(), true)

	Directory d5("dir_f");
	TEST_EQUAL(d5.isValid(), false)
RESULT

CHECK(Directory(const Directory& directory))
	Directory d("dir_a"+PS);
	Directory d1(d);
	TEST_EQUAL(d1 == d, true)
RESULT

CHECK(void clear())
	Directory d("dir_a"+PS);
	d.clear();
	TEST_EQUAL(d.getPath(), "")
RESULT

CHECK(void destroy())
	Directory d1("dir_a"+PS);
	d1.destroy();
RESULT

CHECK(bool set(const String& directory_path, bool set_current = false))
	Directory d;
	d.set("dir_a" + PS + "dir_c" + PS);
	TEST_EQUAL(d.getPath(), test_dir + PS + "dir_a" + PS + "dir_c")
	TEST_EQUAL(d.isCurrent(), false)
	d.set("dir_a" + PS + "dir_c"+ PS, true);
	TEST_EQUAL(d.isCurrent(), true)
	d.setCurrent(test_dir);
RESULT

CHECK(void set(const Directory& directory))
	Directory d("dir_a" + PS);
	Directory d1;
	d1.set(d);
	TEST_EQUAL(d1.getPath(), test_dir + PS + "dir_a")
RESULT

CHECK(Directory& operator = (const Directory& directory))
	Directory d("dir_a" + PS);
	Directory d1;
	d1 = d;
	TEST_EQUAL(d1.getPath(), test_dir + PS + "dir_a")	
RESULT

CHECK(const String& getPath() const)
	STATUS("Creating directory")
	Directory d1("dir_a" + PS);
	STATUS("Done.")
	TEST_EQUAL(d1.getPath(), test_dir + PS + "dir_a")
RESULT

CHECK(void get(Directory& directory) const)
	Directory d("dir_a" + PS);
	Directory d1;
	d.get(d1);
	TEST_EQUAL(d1.getPath(), test_dir + PS + "dir_a")
RESULT

CHECK(bool create(String path, const mode_t& mode = 0777))
	Directory d;
	::chdir(test_dir.c_str());
	Directory d1("dir_a" + PS + "dir_c"+ PS, true);
	TEST_EQUAL(d1.isValid(), true)

	bool ok = cleanup();
	TEST_EQUAL(ok, true)

	bool result = d1.create("test1");
	TEST_EQUAL(result, true)
	TEST_EQUAL(d.setCurrent(test_dir), true)
	Directory d2;
	result = d2.create("dir_a" + PS + "dir_c" + PS + "test2");
	TEST_EQUAL(result, true)
	TEST_EQUAL(d2.isCurrent(), true)
	TEST_EQUAL(d.setCurrent(test_dir), true)
	Directory d3;
	TEST_EQUAL(d3.create(test_dir + PS +"dir_a" + PS + "dir_c" + PS + "test3"), true)

	ok = cleanup();
	TEST_EQUAL(ok, true)
RESULT

CHECK(bool getNextEntry(String& entry))
	bool ok = cleanup();
	TEST_EQUAL(ok, true)	

	Directory d1(".", true);
	TEST_EQUAL(d1.isValid(), true)
	d1.create("test1");
	d1.create("test1" + PS + "test2");

	Directory d2("test1");
	String s;	
	
	// On some systems, the order (or even the presence!) of "." and ".." is not clear.
	// To avoid this, we just skip these entries. These are properties of the
	// Filesystem, not the Directory.
	bool result = d2.getNextEntry(s);
	STATUS("getNextEntry : " << result << " = " << s)
	TEST_EQUAL(result, true)
	bool found_test2 = false;
	while (result == true)
	{
		if ((s != ".") && (s != ".."))
		{
			TEST_EQUAL(s, "test2")
			found_test2 = true;
		}
	  result = d2.getNextEntry(s);
	  STATUS("getNextEntry : " << result << " = " << s)
	}
	TEST_EQUAL(found_test2, true)

	ok = cleanup();
	TEST_EQUAL(ok, true)	
RESULT

CHECK(Size countItems())
{
	Directory d0(test_dir + PS + "dir_a" + PS + "dir_c");
	d0.create("test1");

	Directory d1(test_dir + PS + "dir_a" + PS + "dir_c" + PS + "test1");
	d1.create("test1");
	d1.create("test2");
	d1.create("test3");
	d1.create("test4");
	TEST_EQUAL(d1.countItems(), 4)
}
	bool ok = cleanup();
	TEST_EQUAL(ok, true)
RESULT


CHECK(Size countFiles())
	Directory d(test_dir + PS + "dir_a" + PS + "dir_c");
	TEST_EQUAL(d.countFiles(), 2)
	
	bool result = d.create("test1");
	TEST_EQUAL(result, true)

	Directory d1(test_dir + PS + "dir_a" + PS + "dir_c" + PS + "test1");
	TEST_EQUAL(d1.countFiles(), 0)
RESULT


CHECK(Size countDirectories())
{
	bool ok = cleanup();
	TEST_EQUAL(ok, true)

	Directory d(test_dir + PS + "dir_a" + PS + "dir_c");
	TEST_EQUAL(d.countFiles(), 2)
	
	bool result = d.create("test1");
	TEST_EQUAL(result, true)

	Directory d0(test_dir + PS + "dir_a" + PS + "dir_c");
	Size nr_of_subdirs = 1;
	if (d0.has("CVS"))
	{
		nr_of_subdirs ++;
	}
	TEST_EQUAL(d0.countDirectories(), nr_of_subdirs)
	Directory d1(test_dir + PS + "dir_a" + PS + "dir_c" + PS + "test1");
	TEST_EQUAL(d1.countDirectories(), 0)
}
	bool ok = cleanup();
	TEST_EQUAL(ok, true)
RESULT

CHECK(bool remove(String old_path))
	Directory d;
	d.setCurrent(test_dir);
	Directory d1("dir_a" + PS + "dir_c" + PS, true);
	d1.create("test1");
	TEST_EQUAL(d1.isValid(), true)
	bool result = d1.remove("test1");
	TEST_EQUAL(result, true)
	TEST_EQUAL(d1.remove("xxxx"), false)
RESULT

CHECK(bool rename(String old_path, String new_path))
	bool ok = cleanup();
	TEST_EQUAL(ok, true)

	Directory d;
	d.setCurrent(test_dir);
	Directory d1("dir_a" + PS + "dir_c", true);
	TEST_EQUAL(d1.isValid(), true)
	bool result = d1.create("test1");
	TEST_EQUAL(result, true)	
	result = d1.rename("test1", "test2");
	TEST_EQUAL(result, true)  
	result = d1.rename("test2", "test3");
	TEST_EQUAL(result, true)	  
	TEST_EQUAL(d1.rename("testz", "test"), false)	 
	result = d1.remove("test3");
	TEST_EQUAL(result, true)
	d.setCurrent(test_dir);
RESULT

CHECK(bool renameTo(String new_path))
	Directory d;
	Directory d1("dir_a" + PS + "dir_c");
	d1.remove("test1");
	bool result = d1.create("test1");
	TEST_EQUAL(result, true)
	d.setCurrent(test_dir);
	Directory d2("dir_a" + PS + "dir_c" + PS + "test1", true);
	TEST_EQUAL(d2.isValid(), true)
  result = d2.renameTo("test2");
	TEST_EQUAL(result, true)	  
	result = d1.remove("test2");
	TEST_EQUAL(result, true)
	d.setCurrent(test_dir);
RESULT

CHECK(bool getFirstEntry(String& entry))
	String s;
	Directory d0(test_dir + PS + "dir_a" + PS + "dir_c");

	d0.create("test1");

	TEST_EQUAL(d0.getFirstEntry(s), true)
	TEST_EQUAL((s == ".") || (s == "..") || (s == "test1"), true)
	TEST_EQUAL(d0.getFirstEntry(s), true)
	TEST_EQUAL((s == ".") || (s == "..") || (s == "test1"), true)
	d0.remove("test1");
RESULT

CHECK(bool find(const String& filename, String& filepath))
	Directory d1(test_dir + PS + "dir_a");
	String s;
	bool result = d1.find("dir_c", s);
	TEST_EQUAL(result, true)
	TEST_EQUAL(s, test_dir + PS + "dir_a");
		d1.set(test_dir);
	result = d1.find("a", s);
	TEST_EQUAL(result, true)
	TEST_EQUAL(s, test_dir)
	result = d1.find("f", s);
	TEST_EQUAL(result, false)
RESULT

CHECK(bool has(const String& item))
	Directory d1("dir_a" + PS + "dir_c"+ PS);
	TEST_EQUAL(d1.has("a"), true)
	TEST_EQUAL(d1.has("x"), false)
RESULT

CHECK(bool isCurrent() const)
	Directory d;
	Directory d1("dir_a" + PS + "dir_c"+ PS);
	TEST_EQUAL(d1.isCurrent(), false)
	d1.setCurrent();
	TEST_EQUAL(d1.isCurrent(), true)
	d.setCurrent(test_dir);
RESULT

CHECK(bool isEmpty())
	Directory d1(test_dir + PS + "dir_a" + PS + "dir_c");
	d1.create("xxx");
	Directory d2(test_dir + PS + "dir_a" + PS + "dir_c" + PS + "xxx");
	TEST_EQUAL(d2.isEmpty(), true)
	d1.remove("xxx");
	TEST_EQUAL(d1.isEmpty(), false)
RESULT

CHECK(bool operator == (const Directory& directory) const)
	Directory d("data" + PS);
	Directory d1(d);
	TEST_EQUAL(d == d1, true)
	d1.set("dir_a"+ PS + "dir_c"+ PS);
	TEST_EQUAL(d == d1, false)
RESULT

CHECK(bool operator != (const Directory& directory) const)
	Directory d("data" + PS);
	Directory d1(d);
	TEST_EQUAL(d != d1, false)
	d1.set("dir_a" + PS + "dir_c" + PS);
	TEST_EQUAL(d != d1, true)
RESULT

CHECK(bool isValid() const)
	Directory d("this_directory_should_not_exist");
	TEST_EQUAL(d.isValid(), false)
	Directory a(test_dir + PS);
	TEST_EQUAL(a.isValid(), true)
RESULT

CHECK(bool remove())
	Directory d(test_dir + PS); 
	bool result = d.create("dir_x");
	TEST_EQUAL(result, true)
	d.set(test_dir + PS + "dir_x");
	TEST_EQUAL(d.isValid(), true)
	result = d.remove();
	TEST_EQUAL(result, true)
	TEST_EQUAL(d.isValid(), false)
	d.set("this_directory_should_not_exist");
	result = d.remove();
	TEST_EQUAL(result, false)
RESULT

CHECK(bool setCurrent())
	Directory d(test_dir);
	TEST_EQUAL(d.setCurrent(), true)
	d.set("this_directory_should_not_exist");
	TEST_EQUAL(d.setCurrent(), false)
RESULT

CHECK(String getUserHomeDir())
	TEST_NOT_EQUAL(Directory::getUserHomeDir(), "")
RESULT

CHECK(bool changeToUserHomeDir())
	TEST_NOT_EQUAL(Directory::changeToUserHomeDir(), false)
RESULT

CHECK([Extra]cleanup)
	Directory d(test_dir + PS + "dir_a" + PS + "dir_c");
	d.remove("test1");
	d.remove("test2");
	d.remove("test3");
RESULT

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
