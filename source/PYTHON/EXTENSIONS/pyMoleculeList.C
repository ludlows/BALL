// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// $Id: pyMoleculeList.C,v 1.3 2004/09/22 21:07:05 oliver Exp $
//

#include <BALL/PYTHON/pyMoleculeList.h>
#include <BALL/KERNEL/molecule.h>
#include <BALL/CONCEPT/composite.h>
#include <BALL/KERNEL/fragment.h>
#include <BALL/KERNEL/atomContainer.h>

namespace BALL
{

	PyMoleculeList::PyMoleculeList()
		: List<Molecule*>()
	{
	}

	PyMoleculeList::PyMoleculeList(const PyMoleculeList& new_list)
		: List<Molecule*>(new_list)
  {
	}

	PyMoleculeList::~PyMoleculeList()
		throw()
	{
	}

	PyMoleculeList::PyMoleculeList(const AtomContainer& fragment, bool selected_only )
	{
		set(fragment, selected_only);
	}

	void PyMoleculeList::set(const AtomContainer& fragment, bool selected_only)
	{
		// clear the old contents of the list
		clear();

		// iterate over all molecules
		AtomContainerConstIterator it = fragment.beginAtomContainer();

    for (; +it; ++it)
    {
      const Molecule* molecule = dynamic_cast<const Molecule*>(&*it);
      if ((molecule != 0) && (it->isSelected() || !selected_only))
      {
        // store the molecule pointer in the list
        push_back(const_cast<Molecule*>(molecule));
			}
		}
  }
}
