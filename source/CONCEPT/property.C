// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// $Id: property.C,v 1.31 2003/08/26 09:17:45 oliver Exp $
//

#include <BALL/CONCEPT/property.h>
#include <BALL/CONCEPT/persistenceManager.h>
#include <BALL/CONCEPT/textPersistenceManager.h>

using namespace std;

namespace BALL 
{


	NamedProperty::NamedProperty()
		throw()
		: PersistentObject(),
			type_(NONE),
			name_("")
	{
	}

	NamedProperty::NamedProperty(const NamedProperty& property) 
		throw()
		: PersistentObject(property),
			type_(property.type_),
			name_(property.name_)
	{	
		if (type_ != STRING)
		{
			data_ = property.data_;
		} 
		else 
		{
			data_.s = new string(*(property.data_.s));
		}
	}
 
	void NamedProperty::persistentWrite(PersistenceManager& pm, const char* name) const
		throw(Exception::GeneralException)
	{
		pm.writeObjectHeader(this, name);
			pm.writePrimitive((int)type_, "type_");
			pm.writePrimitive(String(name_), "name_");
			
			switch (type_)
			{
				case	INT:					pm.writePrimitive(data_.i, "data_.i");		break;
				case	FLOAT:				pm.writePrimitive(data_.f, "data_.f");		break;
				case	DOUBLE:				pm.writePrimitive(data_.d, "data_.d");		break;
				case	UNSIGNED_INT:	pm.writePrimitive(data_.ui, "data_.ui"); break;
				case	BOOL:					pm.writePrimitive(data_.b, "data_.b");		break;
				case	OBJECT:				pm.writeObjectPointer(data_.object, "data_.object"); break;
				case	NONE:					break;
				case	STRING:				pm.writePrimitive(String(*data_.s), "data_.s");	break;
				default:
					Log.error() << "cannot write unknown property type: " << (int)type_ << endl;
			}
		pm.writeObjectTrailer(name);
	}
	
	void NamedProperty::persistentRead(PersistenceManager& pm)
		throw(Exception::GeneralException)
	{
		// Clear potentially allocated strings.
		if (type_ == STRING)
		{
			delete data_.s;
		}

		int type;
		pm.readPrimitive(type, "type_");
		type_ = (Type)type;
		String s;
		pm.readPrimitive(s, "name_");
		name_ = s;

		switch (type_)
		{
			case	INT:					pm.readPrimitive(data_.i, "data_.i");		break;
			case	FLOAT:				pm.readPrimitive(data_.f, "data_.f");		break;
			case	DOUBLE:				pm.readPrimitive(data_.d, "data_.d");		break;
			case	UNSIGNED_INT:	pm.readPrimitive(data_.ui, "data_.ui"); break;
			case	BOOL:					pm.readPrimitive(data_.b, "data_.b");		break;
			case	NONE:	break;

			case	STRING:				
				// we have to create a new string
				pm.readPrimitive(s, "data_.s");	
				data_.s = new string(s);
				break;

			case	OBJECT:				
				// the persistence manager will take care of
				// reading the object and will set this pointer afterwards 
				pm.readObjectPointer(data_.object, "data_.object"); 
				break;

			default:
				Log.error() << "Unknown type while reading NamedProperty: " << (int)type_ << endl;
		}
	}

  void PropertyManager::write(PersistenceManager& pm) const
		throw()
  {
		pm.writeStorableObject(bitvector_, "bitvector_");
		Size size = (Size)named_properties_.size();
		pm.writePrimitive(size, "size");
		for (Size i = 0; i < size; i++)
		{
			named_properties_[i].persistentWrite(pm, "");
		}
	}

  bool PropertyManager::read(PersistenceManager& pm)
		throw()
  {
		if (!pm.readStorableObject(bitvector_, "bitvector_"))
		{
			return false;
		}
		
		NamedProperty property("");
		named_properties_.clear();
		Size size = 0;
		pm.readPrimitive(size, "size");
		for (Size i = 0; i < size; i++)
		{
			pm.checkObjectHeader(property, "");
				property.persistentRead(pm);
			pm.checkObjectTrailer("");
			named_properties_.push_back(property);
		}

		return true;
	}
  
	void PropertyManager::set(const PropertyManager& property_manager)
		throw()
	{
		bitvector_ = property_manager.bitvector_;
		named_properties_ = property_manager.named_properties_;
	}

	void PropertyManager::setProperty(const NamedProperty& property)
		throw()
	{
		// search whether the property already exists
		vector<NamedProperty>::iterator it = named_properties_.begin();
		for (; it != named_properties_.end(); ++it)
		{
			if (it->getName() == property.getName())
			{
				// yes, it exists. Erase the old content
				named_properties_.erase(it);
				break;
			}
		}

		// add the new content
		named_properties_.push_back(property);
	}

	void PropertyManager::setProperty(const string& name)
		throw()
	{
		// search whether a property with the same name already exists
		vector<NamedProperty>::iterator it = named_properties_.begin();
		for (; it != named_properties_.end(); ++it)
		{
			if (it->getName() == name)
			{
				// yes, it exists. Erase the old content
				named_properties_.erase(it);
				break;
			}
		}

		// add the new content
		named_properties_.push_back(NamedProperty(name));
	}

	const NamedProperty& PropertyManager::getProperty(const string& name) const
		throw()
	{
		for (Size i = 0; i < named_properties_.size(); ++i)
		{
			if (named_properties_[i].getName() == name)
			{
				return named_properties_[i];
			}
		}

		return RTTI::getDefault<NamedProperty>();
	}
	
	void PropertyManager::clearProperty(const string& name)
		throw()
	{
    vector<NamedProperty>::iterator it = named_properties_.begin();
		for (; it != named_properties_.end(); ++it)
 		{
			if (it->getName() == name)
			{
				named_properties_.erase(it);
				break;
			}
		}
	}

	bool PropertyManager::hasProperty(const string& name) const
		throw()
	{
		for (Size i = 0; i < named_properties_.size(); i++)
		{
			if (named_properties_[i].getName() == name)
			{
				return true;
			}
		}
		
		return false;
	}

	void PropertyManager::dump(ostream& s, Size depth) const
		throw()
	{
		BALL_DUMP_STREAM_PREFIX(s);

		BALL_DUMP_DEPTH(s, depth);
		BALL_DUMP_CLASS_HEADER(s, PropertyManager, this);

		BALL_DUMP_DEPTH(s, depth);
		s << bitvector_ << endl;
		
		for (Size i = 0; i < named_properties_.size(); ++i)
		{
			BALL_DUMP_DEPTH(s, depth);
			s << "(" << named_properties_[i].getName();
			switch (named_properties_[i].getType())
			{
				case NamedProperty::NONE:
					s << "NONE: ";
					break;
				case NamedProperty::BOOL:
					s << "BOOL: " << named_properties_[i].getBool();
					break;
				case NamedProperty::INT:
					s << "INT: " << named_properties_[i].getInt();
					break;
				case NamedProperty::UNSIGNED_INT:
					s << "UNSIGNED_INT: " << named_properties_[i].getUnsignedInt();
					break;
				case NamedProperty::FLOAT:
					s << "FLOAT: " << named_properties_[i].getFloat();
					break;
				case NamedProperty::DOUBLE:
					s << "DOUBLE: " << named_properties_[i].getDouble();
					break;
				case NamedProperty::STRING:
					s << "STRING: " << (char)34 << named_properties_[i].getString().c_str() << (char)34;
					break;
				case NamedProperty::OBJECT:
					s << "OBJECT: " << named_properties_[i].getObject();
			}

			s << ")" << endl;
		}

		BALL_DUMP_STREAM_SUFFIX(s);
	}

#	ifdef BALL_NO_INLINE_FUNCTIONS
#		include <BALL/CONCEPT/property.iC>
#	endif

} // namespace BALL
