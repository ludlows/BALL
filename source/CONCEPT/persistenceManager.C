// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// $Id: persistenceManager.C,v 1.23 2005/02/06 09:45:00 oliver Exp $
//

#include <BALL/CONCEPT/persistenceManager.h>
#include <BALL/KERNEL/system.h>
#include <BALL/KERNEL/protein.h>
#include <BALL/KERNEL/nucleicAcid.h>
#include <BALL/KERNEL/bond.h>

// #define BALL_DEBUG_PERSISTENCE

#ifdef BALL_DEBUG_PERSISTENCE
# define DEBUG(a) Log.info() << a << endl;
#else
# define DEBUG(a)
#endif

using namespace std;

namespace BALL 
{

	PersistenceManager::PersistenceManager()
		throw()
		: create_methods_(),
			object_out_(),
			object_out_needed_(),
			pointer_map_(),
			pointer_list_(),
			object_in_(),
			ostr_(),
			istr_()
	{
		// all kernel classes	
		registerKernelClasses_();
	}

	PersistenceManager::PersistenceManager(const PersistenceManager& pm)
		throw()
		: create_methods_(pm.create_methods_),
			object_out_(pm.object_out_),
			object_out_needed_(pm.object_out_needed_),
			pointer_map_(pm.pointer_map_),
			pointer_list_(pm.pointer_list_),
			object_in_(pm.object_in_),
			ostr_(pm.ostr_),
			istr_(pm.istr_)
	{
	}

	PersistenceManager::PersistenceManager(istream& is)
		throw()
		: create_methods_(),
			object_out_(),
			object_out_needed_(),
			pointer_map_(),
			pointer_list_(),
			object_in_(),
			ostr_(),
			istr_()
	{
		// all kernel classes	
		registerKernelClasses_();
		setIstream(is);
	}

	PersistenceManager::PersistenceManager(ostream& os)
		throw()
		: create_methods_(),
			object_out_(),
			object_out_needed_(),
			pointer_map_(),
			pointer_list_(),
			object_in_(),
			ostr_(),
			istr_()
	{
		// all kernel classes	
		registerKernelClasses_();
		setOstream(os);
	}

	PersistenceManager::PersistenceManager(istream& is, ostream& os)
		throw()
		: create_methods_(),
			object_out_(),
			object_out_needed_(),
			pointer_map_(),
			pointer_list_(),
			object_in_(),
			ostr_(),
			istr_()
	{
		// all kernel classes	
		registerKernelClasses_();
		setIstream(is);
		setOstream(os);
	}

	PersistenceManager::~PersistenceManager()
		throw()
	{
	}

	void PersistenceManager::registerKernelClasses_()
		throw()
	{
		// all kernel classes, their base classes, 
		// and the classes used in kernel classes

		#define REGISTER_CLASS(T) registerClass(RTTI::getStreamName<T>(), T::createDefault);
		REGISTER_CLASS(AtomContainer)
		REGISTER_CLASS(NamedProperty)
		REGISTER_CLASS(Composite)
		REGISTER_CLASS(Atom)
		REGISTER_CLASS(Bond)
		REGISTER_CLASS(Fragment)
		REGISTER_CLASS(System)
		REGISTER_CLASS(Molecule)
		REGISTER_CLASS(PDBAtom)
		REGISTER_CLASS(Residue)
		REGISTER_CLASS(Chain)
		REGISTER_CLASS(Protein)
		REGISTER_CLASS(SecondaryStructure)
		REGISTER_CLASS(NucleicAcid)
		REGISTER_CLASS(Nucleotide)
		#undef REGISTER_CLASS
	}

	PersistenceManager& PersistentObject::operator >> (PersistenceManager& pm) const
		throw(Exception::GeneralException)
	{ 
		DEBUG("PersistenceManager: entering operator >>")
		pm.startOutput();
		persistentWrite(pm, "");
		pm.endOutput();
		return pm;
	}

	void PersistentObject::persistentWrite(PersistenceManager& /* pm */, const char* /* name */) const
		throw(Exception::GeneralException)
	{
	}

	void PersistentObject::persistentRead(PersistenceManager& /* pm */)
		throw(Exception::GeneralException)
	{
	}

	void PersistenceManager::registerClass(String signature, const CreateMethod	m)
		throw()
	{
		create_methods_.insert(signature, m);
	}

	void* PersistenceManager::createObject(String signature) const
		throw()
	{
		DEBUG("PersistenceManager: createObject(" << signature)
		if (create_methods_.has(signature)) 
		{
			CreateMethod method = create_methods_[signature];

			DEBUG("PersistenceManager: created object of type" << signature)
			return method();
		} 
		else 
		{
			return 0;
		}
	}

	Size PersistenceManager::getNumberOfClasses() const
		throw()
	{
		return create_methods_.size();
	}

	void PersistenceManager::setOstream(ostream& s) 
		throw()
	{
		ostr_ = &s;
		object_out_.clear();
		object_out_needed_.clear();
	}

	void PersistenceManager::setIstream(istream& s) 
		throw()
	{
		istr_ = &s;
		pointer_list_.clear();
		pointer_map_.clear();
	}
		
	void PersistenceManager::startOutput()
		throw()
	{
		initializeOutputStream();
		object_out_.clear();
		object_out_needed_.clear();
		writeStreamHeader();
	}

	void PersistenceManager::endOutput()
		throw()
	{
		writeStreamTrailer();
		addNeededObjects_();
		writeStreamTrailer();
		finalizeOutputStream();
	}

	void PersistenceManager::initializeInputStream()
		throw()
	{
	}

	void PersistenceManager::finalizeInputStream()
		throw()
	{
	}

	void PersistenceManager::initializeOutputStream()
		throw()
	{
	}

	void PersistenceManager::finalizeOutputStream()
		throw()
	{
	}

	PersistenceManager& PersistenceManager::operator << (const PersistentObject& object)
		throw()
	{
		object >> *this;
		return *this;
	}

	PersistenceManager& PersistenceManager::operator >> (PersistentObject*& object_ptr)
		throw()
	{
		object_ptr = readObject();
		return *this;
	}
	
	PersistentObject*	PersistenceManager::readObject()
		throw(Exception::GeneralException)
	{
		DEBUG("PersistenceManager: entering readObject")
		if (istr_ == 0)
		{
			return 0;
		}

		PersistentObject*	first_object = 0;
		PersistentObject*	obj = 0;

		pointer_map_.clear();
		pointer_list_.clear();
		object_in_.clear();
		
		String type_name;
		LongSize ptr;

		// prepare the input stream
		initializeInputStream();

		// if an error happened, just exit the loop 
		// to clean up the mess 
		bool error = false;

		// loop while the stream is not empty, 
		// we did not read the END mark and
		// an error did not occur
		while (*istr_ && checkStreamHeader() && !error) 
		{
			// retrieve the first object signature
			getObjectHeader(type_name, ptr);
				
			if (!create_methods_.has(type_name)) 
			{
				// something bad happend - abort everything and clean up!
				Log.error() << "Cannot create object of unregistered class " << type_name << "!" << endl;
				error = true;
				break;
			} 
			
			// Create an instance of type_name 
			CreateMethod	m = create_methods_[type_name];
			obj = (PersistentObject*)m();
			DEBUG("PersistenceManager: created object of type " << type_name << " @ " << (void*)obj)

			// check whether the creation was successful
			if (obj == 0)
			{
				Log.error() << "Could not create object of type " << type_name << "!" << endl;
				error = true;
				break;
			}

			// remember the this pointer of the new object 
			object_in_.push_back(obj);

			// check 
			if (ptr == 0)
			{
				Log.error() << "Read invalid object pointer!" << endl;
				error = true;
				break;
			}
				
			// store the old pointer and the new pointer
			addPointerPair_(ptr, (void*)obj);
				
			// make the new object read itself
			DEBUG("PersistenceManager: calling persistentRead on new object.")
			obj->persistentRead(*this);


			if (!checkObjectTrailer("") || !checkStreamTrailer())
			{
				DEBUG("PersistenceManager: check for object trailer or stream trailer failed")
				error = true;
				break;
			}

			// remember the first object created.
			// this is our return value
			if (first_object == 0)
			{
				first_object = obj;
			}
		}

		// prepare the input stream for closing
		DEBUG("PersistenceManager: calling finalizeInputStream")
		finalizeInputStream();

		if (error)
		{ 
			return 0;
		}

		// update the pointers: replace all old object addresses with the
		// new ones
		if (!updatePointers_())
		{
			return 0;
		}

		// finalize all objects
		ObjectList::iterator object_it = object_in_.begin();
		for (; object_it != object_in_.end(); object_it++)
		{
			const_cast<PersistentObject*>(*object_it)->finalize();
		}

		// return the first object read
		return first_object;
	}

	void PersistenceManager::addPointerPair_(LongSize old_ptr, void* new_ptr)
		throw()
	{
		DEBUG("PersistenceManager: pointer pair (" << hex << old_ptr << "/" << new_ptr << ")")
		pointer_map_.insert(std::make_pair(old_ptr, new_ptr));
	}

	void PersistenceManager::addNeededObjects_() 
		throw(Exception::GeneralException)
	{
		while (object_out_needed_.size() > 0)
		{
			const PersistentObject*	obj = object_out_needed_.back();
			object_out_needed_.pop_back();

			if (!object_out_.has(obj))
			{
				writeStreamHeader();
				(*obj).persistentWrite(*this);
				writeStreamTrailer();
			} 
				
		}
	}

	bool PersistenceManager::updatePointers_()
		throw()
	{
		// assume everything will go smoothly
		bool result = true;

		PointerList::iterator it = pointer_list_.begin();
		for (; it != pointer_list_.end(); ++it) 
		{
			if (pointer_map_.has((*it).second)) 
			{
				// OK. We know the correct value for the pointer
				(*(*it).first) = pointer_map_[(*it).second];
			} 
			else 
			{ 
				Log.error() << "PersistenceManager: size of pointer map: " << pointer_map_.size() << std::endl;
				Log.error() << "PersistenceManager: Could not assign object for pointer to "
																				<< hex << (unsigned int)(*it).second << endl;
				result = false;
			}
		}
		
		return result;
	}

#ifdef BALL_NO_INLINE_FUNCTIONS
#	include <BALL/CONCEPT/persistenceManager.iC>
#endif

} // namespace BALL
