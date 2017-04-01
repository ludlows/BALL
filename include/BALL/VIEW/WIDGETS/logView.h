// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// $Id: logView.h,v 1.11.8.6 2005/11/01 21:47:35 amoll Exp $
//

#ifndef BALL_VIEW_WIDGETS_LOGVIEW_H
#define BALL_VIEW_WIDGETS_LOGVIEW_H

#ifndef BALL_COMMON_H
#	include <BALL/common.h>
#endif

#ifdef BALL_HAS_SSTREAM
# include <sstream>
#else
# include <strstream>
#endif

#ifndef QAPPLICATION_H
#	include <qapplication.h>
#endif

#ifndef QSTRING_H
#	include <qstring.h>
#endif

#ifndef QTEXTEDIT_H
# include <qtextedit.h>
#endif

#ifndef BALL_VIEW_WIDGETS_DOCKWIDGET_H
#	include <BALL/VIEW/WIDGETS/dockWidget.h>
#endif

namespace BALL
{
	namespace VIEW
	{
		class DragLogView
			: public QTextEdit
		{
			Q_OBJECT

			public:

			DragLogView(QWidget* parent);

			public slots:
			virtual void contentsDragEnterEvent(QDragEnterEvent* e);
			virtual void contentsDragLeaveEvent(QDragLeaveEvent* e);
			virtual void contentsDropEvent(QDropEvent* e);
		};

		/** LogView class.
				The class LogView records all messages sent to the  \link BALL::LogStream Log \endlink  object and
				displays them as a text history. The class is derived from 
				<b> NotificationTarget<LogStreamNotifier></b> that provides the connection
				to the  \link BALL::LogStream Log \endlink  object. The class  QTextEdit from the 
				qt - library is responsible for the visualization of the text history.
				Use the class LogView as a widget. There are no initializations necessary.
				Just create this widget as a child widget of your application and it will
				record and show all messages sent to the \link BALL::LogStream Log \endlink object.
				\ingroup ViewWidgets
		*/
		class BALL_VIEW_EXPORT LogView
			: public DockWidget,
			  public LogStreamNotifier
		{
			Q_OBJECT

			public:

			BALL_EMBEDDABLE(LogView,DockWidget)
		
			/**	@name	Constructors
			*/	
			//@{

			/** Default Constructor.
					The contructor connects the own
					<b> stringstream</b> with the  \link BALL::LogStream Log \endlink  object. If a string is written into
					 \link BALL::LogStream Log \endlink  this will be notified and the string will be displayed
					by this logView. 
					\see         BALL::LogStream
			*/
			LogView(QWidget *parent = 0, const char *name = 0)
				throw();

			/** Copy constructor.
				 	Only for Python Interface
					The text of <b> view</b> will be copied into this logView.
			*/
			LogView(const LogView& view)
				throw();

			//@}
			/** @name Destructors */
			//@{

			/** Destructor.
					Calls  clear.
			*/
			virtual ~LogView()
				throw();

			//@}

			/**	Setup the menu entry in "Edit->Clear Logs".
			*/
			virtual void initializeWidget(MainControl& main_control)
				throw();

			/**	Remove menu entries.
			*/
			virtual void finalizeWidget(MainControl& main_control)
				throw();

			protected:

			/** Overridden notify call.
					Will be called by  \link BALL::LogStream Log \endlink  whenever a string is written to it.
					That string will then be displayed.
					\param   source the notification source
					\return  bool returns always <tt>true</tt>
			*/
			void logNotify();

			private:


			QTextEdit* text_edit_;

			bool output_running_;
		};
  	
} } // namespaces

#endif // BALL_VIEW_WIDGETS_LOGVIEW_H
