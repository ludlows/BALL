// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// $Id: pyCompositeList.sip,v 1.2.6.1 2005/11/16 20:01:24 amoll Exp $
//

%MappedType PyCompositeList
{
%TypeHeaderCode
	#include <BALL/PYTHON/pyCompositeList.h>
	#include "sipBALLComposite.h"
  #include <BALL/PYTHON/pyBALLSipHelper.h>
%End

%ConvertFromTypeCode
	PyObject *pl;

	if ((pl = PyList_New(0)) == NULL)
	{
		return NULL;
	}

	// Convert the list.

	for (PyCompositeList::ConstIterator it = sipCpp->begin(); it != sipCpp->end(); ++it)
	{
		Composite& obj = **it;
		PyObject *inst;

		if ((inst = pyMapBALLObjectToSip(obj)) == NULL || PyList_Append(pl,inst) < 0)
		{
			Py_DECREF(pl);
			return NULL;
		}
	}

	return pl;
%End
 
%ConvertToTypeCode
	// Convert a Python list of Composite instances to an CompositeList object on the
	// heap.
	//
	
	if (sipIsErr == NULL)
	{
		return PyList_Check(sipPy);
	}
 
	PyCompositeList *composite_list = new PyCompositeList;
 
	for (int i = 0; i < PyList_GET_SIZE(sipPy); ++i)
	{
		Composite* composite;
 
		composite = reinterpret_cast<Composite*>(sipForceConvertTo_Composite(PyList_GET_ITEM(sipPy,i),sipIsErr));
 
		if (*sipIsErr)
		{
			delete composite_list;
			return 0;
		}
 
		composite_list -> push_back(composite);
	}
 
	*sipCppPtr = composite_list;

	return 1;
%End
};

