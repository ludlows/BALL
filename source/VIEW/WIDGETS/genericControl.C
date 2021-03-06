// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// $Id: genericControl.C,v 1.14.6.3 2005/09/02 14:17:04 amoll Exp $
//

#include <BALL/VIEW/WIDGETS/genericControl.h>
#include <BALL/VIEW/KERNEL/mainControl.h>
#include <BALL/VIEW/KERNEL/message.h>
#include <qlistview.h>

namespace BALL
{
	namespace VIEW
	{


		void MyListView::keyPressEvent(QKeyEvent * e)
		{
			if (e->key() != Key_Delete ||
					parentWidget() == 0 ||
					!RTTI::isKindOf<GenericControl>(*parentWidget())) 
			{
				QListView::keyPressEvent(e);
				return;
			}

			((GenericControl*)parentWidget())->deleteCurrentItems();
		}

		void GenericControl::onContextMenu_(QListViewItem* /*item*/, const QPoint& /*point*/, int /*column*/)
		{}

		GenericControl::GenericControl(QWidget* parent, const char* name)
			throw()
				:	DockWidget(parent, name),
					context_item_(0),
					listview(new MyListView(this))
		{
			// appearance
			listview->setRootIsDecorated(true);
			listview->setSorting(-1);
			listview->setSelectionMode(QListView::Extended);
			setGuest(*listview);	

			connect(listview, SIGNAL(rightButtonPressed(QListViewItem*, const QPoint&, int)), this,
							SLOT(onContextMenu_(QListViewItem*, const QPoint&, int)));
		}


		GenericControl::~GenericControl()
			throw()
		{
			#ifdef BALL_VIEW_DEBUG
				Log.error() << "Destructing object " << this << " of class GenericControl" << endl;
			#endif 

			listview->clear();
		}


		GenericControl::ItemList GenericControl::getSelectedItems()
			throw()
		{
			ItemList selected;
			QListViewItemIterator it(listview);
			for (; it.current(); ++it)
			{
				if (it.current()->isSelected())
				{
					selected.push_back(it.current());
				}
			}

			return selected;
		}

		void GenericControl::deselectOtherControls_()
		{
			notify_(new DeselectControlsMessage);
		}

		void GenericControl::updateSelection() 
		{
			QListViewItemIterator it(listview);
			for (; it.current(); ++it)
			{
				QListViewItem* item = it.current();
				if (item->isSelected())
				{
					deselectOtherControls_();
					return;
				}
			}
		}

		void GenericControl::onNotify(Message *message)
			throw()
		{
			if (!RTTI::isKindOf<DeselectControlsMessage>(*message))
			{
				return;
			}

			listview->selectAll(false);
		}

		void GenericControl::initializeWidget(MainControl& main_control)
			throw()
		{
			main_control.insertDeleteEntry();
			DockWidget::initializeWidget(main_control);
		} 

		void GenericControl::removeItem_(QListViewItem* item, bool update)
			throw()
		{
			while (item->firstChild() != 0)
			{
				removeItem_(item->firstChild(), false);
			}

			delete item;

			if (update) listview->triggerUpdate();
		}

	} // namespace VIEW
} // namespace BALL
