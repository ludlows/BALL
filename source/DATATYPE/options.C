// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// $Id: options.C,v 1.26 2003/08/26 09:17:45 oliver Exp $ 
//

#include <BALL/DATATYPE/options.h>

#include <stdlib.h>
#include <errno.h>
#include <stdio.h>
#include <math.h>

#include <fstream>
#include <list>
#include <algorithm>

using namespace std;
using std::ofstream;
using std::ios;

namespace BALL 
{

	const Size Options::MAX_ENTRY_LENGTH = 1024;

	Options::Options()
	 throw()
		:	StringHashMap<String>(),
			name_("")
	{
	}

	Options::Options(const Options& options)
	 throw()
		:	StringHashMap<String>(options),
			name_(options.name_)
	{
	}

	Options::~Options()
	 throw()
	{
	}

	bool Options::isReal(const String& key) const
	 throw()
	{
		errno = 0;
		char*	endptr;
		String value(get(key));

		// an empty String is no real number
		if (value =="")
		{
			return false;
		}
		
		// try to convert it to a number
		::strtod(value.c_str(), &endptr);

		// return and tell whether it happend to work
		return (errno == 0) && (endptr != value.c_str());
	}

	bool Options::isVector(const String& key) const 
	 throw()
	{
		// if the key does not exist - then the nonexistent value	
		// cannot contain a vector
		if (!has(key))
		{
			return false;
		}

		// try to interpret the string as three double values
		double dummy;
		if (sscanf(get(key).c_str(), "(%lf %lf %lf)", &dummy, &dummy, &dummy) == 3)
		{
			return true;
		}
		
		return false;
	}

	bool Options::isBool(const String& key) const
	 throw()
	{
		String s = get(key);
		if (s == "")
		{
			return false;
		}

		s.toLower();
		return (s.compare("true") == 0 || s.compare("false") == 0);
	}
	
	bool Options::isSet(const String& key) const
	 throw()
	{
		return (StringHashMap<String>::find(key) != StringHashMap<String>::end());
	}
	
	bool Options::isInteger(const String& key) const 
	 throw()
	{
		// if it cannot be read as a floating point number
		// it cannot be an integer
		if (!isReal(key))
		{
			return false;
		}			
		
		// check wheter it is an integer
		long long_value = ::atol(get(key).c_str());
		double double_value   = ::atof(get(key).c_str());

		// check if it is an integer (cutoff is 1e-7)
		if (fabs(double_value - ((double)long_value)) <= 1e-7)
		{
			return true;
		}

		// it is a floating point number, but no integer
		return false;
	}


	double Options::getReal(const String& key) const 
	 throw()
	{
		if (!has(key))
		{
			return 0.0;
		}
		double value;
		errno = 0;
		value = ::atof((*find(key)).second.c_str());
		
		if (errno == 0)
		{
			return value;
		} 
		else 
		{
			return 0.0;
		}
	}

	Vector3	Options::getVector(const String& key) const 
	 throw()
	{
		Vector3	h(0,0,0);
		
		if (!has(key))
		{
			return h;
		}
		
		// use temporary variables of double precision
		// this avoids trouble if the definition of 
		// Vector3 is changed between TVector3<float> and TVector3<double>
		double x, y, z;
		
		// try to interpret the option as a vector
		sscanf(get(key).c_str(), "(%lf %lf %lf)", &x, &y, &z);
		h.x = x; 
		h.y = y;
		h.z = z;

		return h;
	}


	bool Options::getBool(const String& key) const 
	 throw()
	{
		if (!has(key))
		{
			return false;
		}

		ConstIterator it = find(key);
		if ((it != end()) && (it->second == "true"))
		{
			return true;
		}
		else 
		{
			return false;
		}
	}

	long Options::getInteger(const String& key) const 
	 throw()
	{
		if (!has(key))
		{
			return 0;
		}

		long value;
		errno = 0;
		ConstIterator it = find(key);
		if (it == end())
		{
			return 0;
		}
		value = atol((*it).second.c_str());
		
		if (errno == 0)
		{
			return value;
		} 
		else 
		{
			errno = 0;
			return 0;
		}
	}

	void Options::set(const String& key, const String& value)
	 throw()
	{
		(*this)[key] = value;
	}

	void Options::setInteger(const String& key, const long value)
	 throw()
	{
		static char buffer[MAX_ENTRY_LENGTH + 1];
		sprintf(buffer, "%ld", value);			
		set(key, &(buffer[0]));
	}

	void Options::setReal(const String& key, const double value)
	 throw()
	{
		char buffer[MAX_ENTRY_LENGTH + 1];
		sprintf(buffer, "%f", value);
		set(key, buffer);
	}

	void Options::setVector(const String& key, const Vector3& value)
	 throw()
	{
		char buffer[MAX_ENTRY_LENGTH + 1];
		sprintf(buffer, "(%f %f %f)", value.x, value.y, value.z);
		set(key, buffer);
	}

	void Options::setBool(const String& key, const bool value)
	 throw()
	{
		if (value)
		{
			set(key, "true");
		}
		else
		{
			set(key, "false");
		}
	}

	String Options::setDefault(const String& key, const String& value)
	 throw()
	{
		if (!has(key))
		{
			set(key, value);
			return value;
		} 
		else 
		{
			return get(key);
		}
	}
		
	double Options::setDefaultReal(const String& key, const double value)
	 throw()
	{
		if (!has(key) || !isReal(key))
		{
			setReal(key, value);
			return value;
		} 
		else 
		{
			return getReal(key);
		}
	}

	bool Options::setDefaultBool(const String& key, const bool value)
	 throw()
	{
		if (!has(key) || !isBool(key))
		{
			setBool(key, value);
			return value;
		} 
		else 
		{
			return getBool(key);
		}
	}

	long Options::setDefaultInteger(const String& key, const long value)
	 throw()
	{
		if (!has(key) || !isInteger(key))
		{
			setInteger(key, value);
			return value;
		}
		else
		{
			return getInteger(key);
		}
	}

	void Options::setName(const String& name)
	 throw()
	{
		name_ = name;
	}

	const String& Options::getName() const
	 throw()
	{
		return name_;
	}

	String Options::get(const String& key) const
	 throw()
	{
		if (!has(key))
		{
			return "";
		}
		ConstIterator it = find(key);

		if (it == end())
		{
			return "";
		}
		else
		{
			return (*it).second;
		}
	}

	bool Options::readOptionFile(const String& filename)
	 throw()
	{
		ifstream	infile;
		infile.open(filename.c_str(), ios::in);
		if (!infile)
		{
			return false;
		}

		char		buffer[MAX_ENTRY_LENGTH + 1];
		String	s, key;
		while (infile.getline(buffer, MAX_ENTRY_LENGTH))
		{
			if ((buffer[0] != '#') && (buffer[0] != '!') && (buffer[0] != ';')) 
			{
				s = buffer;
				key = s.getField(0, " ");
				s = s.after(" ");
				set(key, s);
			}
		}
					
		infile.close();
		return true;
	}
		
	bool Options::writeOptionFile(const String& filename) const
	 throw()
	{
		std::list<String>		entry_list;
		String							entry;

		std::ofstream stream(filename.c_str());//, File::OUT);
		if (!stream.is_open())
		{
			return false;
		}

		stream << "![OptionsTable: " << getName() << " (" << size() << " entries)]" << endl;

		StringHashMap<String>::ConstIterator	it(begin());
		for(; !(it == end()); ++it)
		{
			entry = (*it).first + ' ' + (*it).second;
			entry_list.push_back(entry);
		}

		entry_list.sort();

		std::list<String>::iterator	list_it = entry_list.begin();
		for (; list_it != entry_list.end(); ++list_it) 
		{
	  	stream << *list_it << endl;
		}

		stream << "!-----------------------------------" << endl;
		
		stream.close();
		entry_list.clear();
		return true;
	}

	void Options::dump(ostream& stream, Size /* depth */) const
	 throw()
	{
		std::list<String>		entry_list;
		String							entry;

		stream << "[OptionsTable: " << getName() << " (" << size() << " entries)]" << endl;

		StringHashMap<String>::ConstIterator	it(begin());
		for(; !(it == end()); ++it)
		{
			entry = (*it).first + ' ' + (*it).second;
			entry_list.push_back(entry);
		}

		entry_list.sort();

		std::list<String>::iterator	list_it = entry_list.begin();
		for (; list_it != entry_list.end(); ++list_it) 
		{
	  	stream << *list_it << endl;
		}

		stream << "-----------------------------------" << endl;

		entry_list.clear();
	}

	const Options& Options::operator = (const Options& options)
	 throw()
	{
		StringHashMap<String>::operator = (options);
		name_ = options.name_;

		return *this;
	}

	bool Options::operator == (const Options& option) const 
		throw()
	{
		if (this->name_ != option.name_)
		{
			return false;
		}
		
		return StringHashMap<String>::operator == (option);
	}

	bool Options::operator != (const Options& option) const 
		throw()
	{
		return !(*this == option);
	}

	void Options::clear() throw()
	{
		name_ = "";
		StringHashMap<String>::clear();
	}


} // namespace BALL
