// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// $Id: dockWidget.C,v 1.24.6.6 2005/11/17 14:28:00 amoll Exp $

#include <BALL/VIEW/WIDGETS/dockWidget.h>
#include <BALL/VIEW/KERNEL/mainControl.h>
#include <BALL/VIEW/KERNEL/message.h>

#include <qmenubar.h>
#include <qlabel.h>
#include <qdragobject.h>

namespace BALL
{
	namespace VIEW
	{

		DockWidget::DockWidget()
			:	QDockWindow(QDockWindow::InDock, 0),
				ModularWidget("unnamed DockWidget"),
				guest_(0)
		{
			setAcceptDrops(true);
		}
			
		DockWidget::DockWidget(const DockWidget& dw)
			:	QDockWindow(QDockWindow::InDock, 0),
				ModularWidget(dw.name()),
				guest_(0)
		{
			setAcceptDrops(true);
		}
			
		
		DockWidget::DockWidget(QWidget* parent, const char* name)
			: QDockWindow(QDockWindow::InDock, parent),
				ModularWidget(name),
				guest_(0)
		{
			layout_ = new QVBoxLayout();
			caption_label_ = new QLabel(this, "caption_label");
			caption_label_->resize(90, 12);
			caption_label_->setSizePolicy( QSizePolicy(QSizePolicy::Expanding, QSizePolicy::Fixed, 0, 0, false));
			caption_label_->setPaletteBackgroundColor( QColor( 255, 255, 127 ) );
			QFont caption_label_font(caption_label_->font());
			caption_label_font.setFamily( "Helvetica" );
			caption_label_font.setPointSize( 11 );
			caption_label_->setFont( caption_label_font ); 
			caption_label_->setFrameShape(QLabel::NoFrame);
			caption_label_->setAlignment(QLabel::AlignCenter);
			layout_->addWidget(caption_label_);
			boxLayout()->addLayout(layout_);

			setOrientation(Qt::Vertical);

			if (name != 0) 
			{ 
				setCaption(name);
				setName(name);
				caption_label_->setText(name);
			}
			else 
			{
				setName("DockWidget");
			}

			setAcceptDrops(true);
		}

		void DockWidget::setGuest(QWidget& guest)
		{
			QPoint p;
			guest.reparent(this, p, true);
			guest.resize(120,100);
			guest.setSizePolicy( QSizePolicy(QSizePolicy::Expanding, QSizePolicy::Expanding, 0, 0, false));
			layout_->addWidget(&guest);
			guest_ = &guest;
			setMinimumSize(20, 20);
			setCloseMode(QDockWindow::Always);
			setResizeEnabled(true);
		}

		void DockWidget::initializeWidget(MainControl&)
			throw()
		{
			window_menu_entry_id_ = 
				insertMenuEntry(MainControl::WINDOWS, getIdentifier(), this, SLOT(switchShowWidget()));
			menuBar()->setItemChecked(window_menu_entry_id_, true);
			connect(this, SIGNAL(visibilityChanged(bool)), this, SLOT(setWindowsMenuEntry(bool)));
		}

		void DockWidget::writePreferences(INIFile& /* inifile */)
			throw()
		{
			// prevent call of ModularWidget::writePreferences()
		}

		void DockWidget::fetchPreferences(INIFile & inifile)
			throw()
		{
			// if the INIFile does not have the information to restore the state of the dockwidgets,
			// make only the default widgets visible
			if (!inifile.hasEntry("WINDOWS", "Main::dockwidgets") &&
				  !default_visible_)
			{
				switchShowWidget();
			}

			if (!BALL_VIEW_DOCKWINDOWS_SHOW_LABELS)
			{
				caption_label_->hide();
			}
		}

		void DockWidget::switchShowWidget()
			throw()
		{
			if (window_menu_entry_id_ == -1) return;

			if (!getMainControl()) 
			{
				BALLVIEW_DEBUG
				return;
			}

			QMenuBar* menu = getMainControl()->menuBar();
			if (menu->isItemChecked(window_menu_entry_id_))
			{
				hide();
				menu->setItemChecked(window_menu_entry_id_, false);
			}
			else
			{
				show();
				menu->setItemChecked(window_menu_entry_id_, true);
			}
		}

		void DockWidget::setWindowsMenuEntry(bool state) 
		{
			QMenuBar* menu = getMainControl()->menuBar();
			if (menu->isItemChecked(window_menu_entry_id_) != isVisible())
			{
				menu->setItemChecked(window_menu_entry_id_, state);
			}
		}

		void DockWidget::applyPreferences()
			throw()
		{
			if (!BALL_VIEW_DOCKWINDOWS_SHOW_LABELS) caption_label_->hide();
			else caption_label_->show();
		}

		void DockWidget::dropEvent(QDropEvent* e)
		{
			VIEW::processDropEvent(e);
		}

		void DockWidget::dragEnterEvent(QDragEnterEvent* event)
		{
			event->accept(QTextDrag::canDecode(event));
		}

		void DockWidget::setVisible(bool state)
		{
			if (state)
			{
				show();
			}
			else
			{
				hide();
			}
		}

		void DockWidget::registerWidgetForHelpSystem(const QWidget* widget, const String& url)
		{
			ModularWidget::registerWidgetForHelpSystem(widget, url);

			// also register Windows menu entry!
			RegisterHelpSystemMessage* msg = new RegisterHelpSystemMessage();
			msg->setMenuEntry(window_menu_entry_id_);
			msg->setURL(url);
			notify_(msg);
		}


	} // namespace VIEW 
} // namespace BALL
