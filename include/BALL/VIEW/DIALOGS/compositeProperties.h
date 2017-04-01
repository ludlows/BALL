// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// $Id: compositeProperties.h,v 1.3.6.1 2005/09/01 22:17:44 amoll Exp $
//

#ifndef BALL_VIEW_DIALOGS_COMPOSITEPROPERTIES_H
#define BALL_VIEW_DIALOGS_COMPOSITEPROPERTIES_H

#include <BALL/VIEW/UIC/compositePropertiesData.h>

#ifndef BALL_CONCEPT_COMPOSITE
# include <BALL/CONCEPT/composite.h>
#endif

namespace BALL
{
	namespace VIEW
	{

		/** Dialog for showing and changing the properties of an Composite, e.g. an Atom or a Molecule.
				\ingroup  ViewDialogs
		*/
		class BALL_VIEW_EXPORT CompositeProperties : public CompositePropertiesData
		{ 
			Q_OBJECT

		public:
			CompositeProperties(Composite* composite,  QWidget* parent = 0, const char* name = 0, 
													bool modal = FALSE, WFlags fl = 0 );
			~CompositeProperties();

		public slots:
			void accept();

		private:
			String getString_(float data) const
				throw();

			Composite* composite_;
		};

	} // namespace VIEW
} // namespace BALL

#endif // BALL_VIEW_DIALOGS_COMPOSITEPROPERTIES_H
