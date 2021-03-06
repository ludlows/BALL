// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// $Id: pyCompositeList.C,v 1.3.6.1 2005/06/15 11:19:02 amoll Exp $
//

#include <BALL/PYTHON/pyCompositeList.h>
#include <BALL/CONCEPT/composite.h>

namespace BALL
{

	PyCompositeList::PyCompositeList()
		: List<Composite*>()
	{
	}

	PyCompositeList::PyCompositeList(const PyCompositeList& composite_list)
		: List<Composite*>(composite_list)
  {
	}

	PyCompositeList::PyCompositeList(const List<Composite*>& composite_list)
		: List<Composite*>(composite_list)
  {
	}

	PyCompositeList::PyCompositeList(const List<const Composite*>& composite_list)
		: List<Composite*>()
  {
		List<const Composite*>::const_iterator it(composite_list.begin());
		for (; it != composite_list.end(); ++it)
		{
			push_back(const_cast<Composite*>(*it));
		}
	}

	PyCompositeList::PyCompositeList(const HashSet<Composite*>& composite_set)
		: List<Composite*>()
  {
		std::copy(composite_set.begin(), composite_set.end(), std::inserter(*this, begin()));
	}

	PyCompositeList::PyCompositeList(const HashSet<const Composite*>& composite_set)
		: List<Composite*>()
  {
		HashSet<const Composite*>::ConstIterator it(composite_set.begin());
		for (; +it; ++it)
		{
			push_back(const_cast<Composite*>(*it));
		}
	}

	PyCompositeList::~PyCompositeList()
		throw()
	{
	}

}
