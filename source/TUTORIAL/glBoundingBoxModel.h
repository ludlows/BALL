// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// $Id: glBoundingBoxModel.h,v 1.4 2004/04/16 11:49:30 amoll Exp $
// BALLView tutorial example
// ------------------------
// create a processor that computes the bounding box of a molecular structure

// BALL includes
#include <BALL/MATHS/simpleBox3.h>

// BALLVIEW includes
#include <BALL/VIEW/PRIMITIVES/simpleBox.h>
#include <BALL/VIEW/MODELS/modelProcessor.h>

namespace BALL
{

	namespace VIEW
	{

		class GLBoundingBoxModel
			: public ModelProcessor
		{
			public:

			// constructors and destructor
			GLBoundingBoxModel() throw();
			virtual ~GLBoundingBoxModel()	throw();

			// accessors
			void setColor(const ColorRGBA& color) throw();

			// processor related methods
			virtual bool start() throw();
			virtual bool finish()	throw();
			virtual Processor::Result operator() (Composite& composite)	throw();

			protected:
			ColorRGBA color_;
			bool new_start_;
			Composite* start_composite_;
			SimpleBox3 bbox_;
		};
	
	}

}
