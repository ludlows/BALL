// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
#ifndef BALL_VIEW_DIALOGS_BONDPROPERTIES_H
#define BALL_VIEW_DIALOGS_BONDPROPERTIES_H

#include <BALL/VIEW/UIC/bondPropertiesData.h>

#ifndef BALL_KERNEL_ATOM_H
# include <BALL/KERNEL/atom.h>
#endif

namespace BALL
{
	namespace VIEW
	{

		/** Dialog for showing and changing the properties of the bonds of an atom
				\ingroup  ViewDialogs
		*/
		class BALL_VIEW_EXPORT BondProperties 
			: public BondPropertiesData
		{ 
			Q_OBJECT

		public:
			BondProperties(Atom* atom,  QWidget* parent = 0, const char* name = 0, bool modal = FALSE, WFlags fl = 0 );
			~BondProperties();

		public slots:
			void bondSelected();
			void focusAtom();
			void focusPartner();

		private:
			Atom* atom_;
			QWidget* parent_;
		};

	} // namespace VIEW

} // namespace BALL

#endif // BALL_VIEW_DIALOGS_BONDPROPERTIES_H
