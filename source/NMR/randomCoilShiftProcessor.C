// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// $Id: randomCoilShiftProcessor.C,v 1.6 2002/02/27 12:23:55 sturm Exp $

#include<BALL/NMR/randomCoilShiftProcessor.h>
#include<BALL/FORMAT/parameterSection.h>

namespace BALL
{

	const char* RandomCoilShiftProcessor::PROPERTY__RANDOM_COIL_SHIFT = "RandomCoilShift";

	RandomCoilShiftProcessor::RandomCoilShiftProcessor()
		throw()
	{
	}

	RandomCoilShiftProcessor::RandomCoilShiftProcessor(const RandomCoilShiftProcessor& processor)
		throw()
		:	ShiftModule(processor),
			shift_map_(processor.shift_map_)
	{
	}

	RandomCoilShiftProcessor::~RandomCoilShiftProcessor()
		throw()
	{
	}

	void RandomCoilShiftProcessor::init()
		throw()
	{
		// we assume the worst case: init fails -> valid_ = false
		valid_ = false;

		// make sure we have correct parameters
		if (parameters_ == 0)
		{
			return;
		}

		// extract the data from the correct section
		ParameterSection parameter_section;
		parameter_section.extractSection(*parameters_, "RandomCoilShifts");

		// make sure the section contains the required columns
		if (!parameter_section.hasVariable("shift"))
		{
			return;
		}

		// extract the random coil shifts from the parameter section
		Position shift_column = parameter_section.getColumnIndex("shift");
		for (Position i = 0; i < parameter_section.getNumberOfKeys(); i++)
		{
			shift_map_.insert(parameter_section.getKey(i), 
												parameter_section.getValue(i, shift_column).toFloat());
		}
				
		// initialization complete - mark the module as valid
		valid_ = true;

		return;
	}

	Processor::Result RandomCoilShiftProcessor::operator () (Composite& composite)
		throw()
	{
		Atom* atom_ptr = dynamic_cast<Atom*>(&composite);
		if (atom_ptr == 0)
		{
			return Processor::CONTINUE;
		}

		String full_name = atom_ptr->getFullName();
		full_name.substitute(":", " ");
		if (!shift_map_.has(full_name))
		{
			full_name = atom_ptr->getFullName(Atom::NO_VARIANT_EXTENSIONS);	
			full_name.substitute(":", " ");
			if (!shift_map_.has(full_name))
			{
				full_name = "* " + atom_ptr->getName();
				if (!shift_map_.has(full_name))
				{
					full_name = "";
				}
			}
		}

		if (full_name != "")
		{
			// retrieve the random coil shift from the hash map
			const float& delta_RC = shift_map_[full_name];

			// add the random coil shift to the total shift value
			float delta = atom_ptr->getProperty(ShiftModule::PROPERTY__SHIFT).getFloat();
			delta += delta_RC;
			atom_ptr->setProperty(ShiftModule::PROPERTY__SHIFT, delta);
			
			// store the random coil shift in the random coil shift property
			atom_ptr->setProperty(RandomCoilShiftProcessor::PROPERTY__RANDOM_COIL_SHIFT, delta_RC);
		}
		
		return Processor::CONTINUE;
	}

}//namespace BALL
