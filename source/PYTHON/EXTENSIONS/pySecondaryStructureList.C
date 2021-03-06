// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// $Id: pySecondaryStructureList.C,v 1.2 2002/02/27 12:24:01 sturm Exp $

#include <BALL/PYTHON/pySecondaryStructureList.h>
#include <BALL/KERNEL/secondaryStructure.h>
#include <BALL/KERNEL/atomContainer.h>

namespace BALL
{

	PySecondaryStructureList::PySecondaryStructureList()
		: List<SecondaryStructure*>()
	{
	}

	PySecondaryStructureList::PySecondaryStructureList(const PySecondaryStructureList& new_list)
		: List<SecondaryStructure*>(new_list)
  {
	}

	PySecondaryStructureList::~PySecondaryStructureList()
		throw()
	{
	}

	PySecondaryStructureList::PySecondaryStructureList(const AtomContainer& fragment, bool selected_only )
	{
		set(fragment, selected_only);
	}

	void PySecondaryStructureList::set(const AtomContainer& fragment, bool selected_only)
	{
		// clear the old contents of the list
		clear();

		// iterate over all secondary structures
		AtomContainerConstIterator it = fragment.beginAtomContainer();

    for (; +it; ++it)
    {
      const SecondaryStructure* sec_struct = dynamic_cast<const SecondaryStructure*>(&*it);
      if ((sec_struct != 0) && (it->isSelected() || !selected_only))
      {
        // store the sec_struct pointer in the list
        push_back(const_cast<SecondaryStructure*>(sec_struct));
			}
		}
	}
}
