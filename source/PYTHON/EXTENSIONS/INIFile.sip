// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// $Id: INIFile.sip,v 1.4 2003/11/01 09:27:36 oliver Exp $
//


class INIFile
{
%TypeHeaderCode
  #include <BALL/FORMAT/INIFile.h>
  using namespace BALL;
%End
	public:
	enum
	{
					MAX_LINE_LENGTH = 1024
	};


	static const String UNDEFINED;
	static const String HEADER;

	INIFile();
	INIFile(const String& /filename/);
	~INIFile();
	void clear();

	bool read();

	bool write();
	const String& getFilename() const;
	void setFilename(const String& /filename/);

	bool isValid() const;

	//LineIterator getLine(Size /line_number/);

	//bool setLine(LineIterator line_it, const String& line);
	//bool deleteLine(LineIterator line_it);

	//bool insertLine(LineIterator line_it, const String& line);

	//bool appendLine(const String& section_name, const String& line);
	Size getNumberOfLines() const;

	bool hasSection(const String& /section_name/) const;
	//Section_iterator getSection(const String& section_name);
	//Section_iterator getSection(Position pos);
	Size getNumberOfSections() const;

	//LineIterator getSectionFirstLine(const String& /section_name/);

	//LineIterator getSectionLastLine(const String& /section_name/);

	Size getSectionLength(const String& /section_name/) const;
	bool deleteSection(const String& /section/);
	bool appendSection(const String& /section/);

	bool hasEntry(const String& /section/, const String& /key/) const;

	String getValue(const String& /section/, const String& /key/) const;

	bool setValue(const String& /section/, const String& /key/, const String& /value/);
	bool operator == (const INIFile& /inifile/) const;
	// bool isValid(const LineIterator& it) const;
	// bool isValid(const Section_iterator& it) const;
	// bool apply(UnaryProcessor<LineIterator>& processor);
	void setDuplicateKeyCheck(bool /mode/);
};
