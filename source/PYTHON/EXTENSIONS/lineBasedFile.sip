// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// $Id: lineBasedFile.sip,v 1.11 2003/11/12 18:20:28 oliver Exp $
//

class LineBasedFile
	:  File
{
%TypeHeaderCode
	#include <BALL/FORMAT/lineBasedFile.h>	
	using namespace BALL;
	using namespace BALL::Exception;
%End
				public:

				LineBasedFile()
								throw();
				LineBasedFile(const String& /filename/, OpenMode /open_mode/ = std::ios::in)
								throw(FileNotFound);
				LineBasedFile(const LineBasedFile& /f/)
								throw(FileNotFound);
				void clear()
								throw();

				bool operator == (const LineBasedFile& /f/) throw();
				bool operator != (const LineBasedFile& /f/) throw();

				//const LineBasedFile& operator = (const LineBasedFile& /file/)
				//				throw();
				Position getLineNumber()
								const throw();
				const String& getLine()
								const throw();
				String& getLine()
								throw();
				bool readLine()
								throw(ParseError);
				bool skipLines(Size /number/ = 1)
								throw(ParseError);

				bool search(const String& /text/, bool /return_to_point/ = false)
								throw(ParseError);

				bool search(const String& /text/, const String& /stop/, 
										bool /return_to_point/ = false)
								throw(ParseError);
				bool gotoLine(Position /line_number/)
								throw(ParseError);
				void rewind()
								throw(ParseError);

				void test(const char* /file/, int /line/, bool /condition/, const String& /msg/)
								const throw(ParseError);

				//String getField(Index /pos/ = 0, const String& /quotes/ = "",
				//								const String& /delimiters/ = String::CHARACTER_CLASS__WHITESPACE)
				//				const throw(IndexUnderflow);

				bool startsWith(const String& /text/)
								const throw();

				bool has(const String& /text/)
								const throw();

				//Index switchString(const std::vector<String>& /data/)
				//				const throw();

				// bool parseColumnFormat(const char* /format/, Position /index/, Size /length/, void* /arg/);
};
