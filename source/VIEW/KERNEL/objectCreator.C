// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// $Id: objectCreator.C,v 1.2 2003/08/26 15:26:32 amoll Exp $

#include <BALL/VIEW/KERNEL/objectCreator.h>

using namespace std;

namespace BALL
{
  
	namespace VIEW
	{
  
		ObjectCreator::ObjectCreator()
				throw()
			:
			init_(false),
			pm_()
		{
		}

		ObjectCreator::~ObjectCreator()
				throw()
		{
			#ifdef BALL_VIEW_DEBUG
				cout << "Destructing object " << (void *)this 
					<< " of class " << RTTI::getName<ObjectCreator>() << endl;
			#endif 
		}

		void ObjectCreator::clear()
				throw()
		{
		}


 	  void ObjectCreator::initPersistenceManager(TextPersistenceManager & /* pm */)
				throw()
    {
    }

	  Composite *ObjectCreator::convertObject(PersistentObject & /* po */)
				throw()
    {
			return (Composite *)0;
    }

	  Composite *ObjectCreator::operator() (IOStreamSocket &iostream_socket)
				throw()
    {
			// initialize the PersistenceManager only one times
			if (init_ == false)
			{
				init_ = true;

				initPersistenceManager(pm_);
			}

			// read persistent object from stream
			pm_.setIstream(iostream_socket);
			PersistentObject*	po = pm_.readObject();
			
			// convert the object
			return convertObject(*po);
    }

	} // namespace VIEW

} // namespace BALL
