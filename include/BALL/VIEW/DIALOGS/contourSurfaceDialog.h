// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
#ifndef BALL_VIEW_DIALOGS_CONTOURSURFACEDIALOG_H
#define BALL_VIEW_DIALOGS_CONTOURSURFACEDIALOG_H

#include <BALL/VIEW/UIC/contourSurfaceDialogData.h>

#ifndef BALL_DATATYPE_REGULARDATA3D_H
 #include <BALL/DATATYPE/regularData3D.h>
#endif

#ifndef BALL_VIEW_DATATYPE_COLORRGBA_H
 #include <BALL/VIEW/DATATYPE/colorRGBA.h>
#endif

namespace BALL
{
 namespace VIEW
 {
	 class DatasetControl;

/** Dialog for creating contour surfaces from RegularData3D
    \ingroup  ViewDialogs
 */
class BALL_VIEW_EXPORT ContourSurfaceDialog 
	: public ContourSurfaceDialogData
{ 
    Q_OBJECT

	public:
    ContourSurfaceDialog( QWidget* parent = 0, const char* name = 0);
    ~ContourSurfaceDialog();

    RegularData3D* getGrid();
    double getThreshold() const;
		void setDatasetControl(DatasetControl* control) {control_ = control;}
		ColorRGBA getColor();

	public slots:
		virtual bool exec();

	protected slots: 
		virtual void valuesChanged();
		virtual void chooseColor();
	
	private:
		RegularData3D* grid_;
		DatasetControl* control_;
};

} } // namespaces
#endif // CONTOURSURFACEDIALOG_H
