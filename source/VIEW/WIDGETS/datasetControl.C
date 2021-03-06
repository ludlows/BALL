// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// $Id: datasetControl.C,v 1.37.2.18 2005/12/18 23:20:27 amoll Exp $
//

#include <BALL/VIEW/WIDGETS/datasetControl.h>
#include <BALL/VIEW/KERNEL/mainControl.h>
#include <BALL/VIEW/KERNEL/message.h>

#include <BALL/VIEW/DIALOGS/snapShotVisualisation.h>
#include <BALL/VIEW/DIALOGS/contourSurfaceDialog.h>

#include <BALL/VIEW/WIDGETS/regularData1DWidget.h>
#include <BALL/VIEW/WIDGETS/regularData2DWidget.h>

#include <BALL/VIEW/PRIMITIVES/mesh.h>
#include <BALL/VIEW/PRIMITIVES/line.h>

#include <BALL/FORMAT/DCDFile.h>
#include <BALL/MOLMEC/COMMON/snapShotManager.h>
#include <BALL/DATATYPE/contourSurface.h>

#include <qfiledialog.h>
#include <qlistview.h>

namespace BALL
{
	namespace VIEW
	{

		DatasetControl::DatasetControl(QWidget* parent, const char* name)
			throw()
			:	GenericControl(parent, name),
				dialog_(0),
				surface_dialog_(0)
		{
		#ifdef BALL_VIEW_DEBUG
			Log.error() << "new DatasetControl " << this << std::endl;
		#endif
			listview->addColumn("Name");
			listview->addColumn("from");
			listview->addColumn("Type");
			listview->setColumnWidth(0, 120);
			listview->setColumnWidth(1, 60);
			listview->setColumnWidth(2, 60);
			default_visible_ = false;
			connect(listview, SIGNAL(selectionChanged()), this, SLOT(updateSelection()));
			registerWidget(this);
		}


		DatasetControl::~DatasetControl()
			throw()
		{
			#ifdef BALL_VIEW_DEBUG
				Log.error() << "Destructing object " << this << " of class DatasetControl" << std::endl;
			#endif 

			if (dialog_ != 0) delete dialog_;

			while (listview->firstChild() != 0)
			{
				deleteItem_(*listview->firstChild());
			}
		}


		void DatasetControl::initializeWidget(MainControl& main_control)
			throw()
		{
			open_trajectory_id_ = insertMenuEntry(MainControl::FILE_OPEN, "Trajectory", 
														this, SLOT(addTrajectory()));
			setMenuHint("Open a trajectory file (1 System has to be selected)");

			insertMenuEntry(MainControl::FILE_OPEN, "1D Grid", this, SLOT(add1DGrid()));
			setMenuHint("Open a 1D data grid");

			insertMenuEntry(MainControl::FILE_OPEN, "2D Grid", this, SLOT(add2DGrid()));
			setMenuHint("Open a 2D data grid");

			insertMenuEntry(MainControl::FILE_OPEN, "3D Grid", this, SLOT(add3DGrid()));
			setMenuHint("Open a 3D data grid");

			menu_cs_ = insertMenuEntry(MainControl::TOOLS, "Contour S&urface", this,  
																							SLOT(computeIsoContourSurface()), CTRL+Key_U);
			setMenuHint("Calculate an isocontour surface from a 3D grid. The grid has to be loaded in the DatasetControl.");
			setMenuHelp("datasetControl.html#isocontour_surfaces");

			GenericControl::initializeWidget(main_control);

			registerWidgetForHelpSystem(this, "datasetControl.html");
		}


		void DatasetControl::checkMenu(MainControl& main_control)
			throw()
		{
			menuBar()->setItemEnabled(open_trajectory_id_, main_control.getSelectedSystem());
			if (getSelectedItems().size() > 0) main_control.setDeleteEntryEnabled(true);

			menuBar()->setItemEnabled(menu_cs_, 
												!getMainControl()->compositesAreLocked() && item_to_grid3_.size() > 0);
		}


		void DatasetControl::addTrajectory()
			throw()
		{
			if (!getMainControl()->getSelectedSystem()) return;

			QString file = QFileDialog::getOpenFileName(
													getWorkingDir().c_str(),
													"DCD files(*.dcd)",
													this,
													"Trajectory File Dialog",
													"Select a DCD file" );

			if (file == QString::null) return;

			addTrajectory(file.ascii());
		}

		void DatasetControl::addTrajectory(const String& filename)
		{
			if (getMainControl()->getSelectedSystem() == 0) return;

			// construct a name for the system(the filename without the dir path)
			DCDFile* dcd = new DCDFile(filename, std::ios::in);
			insertTrajectory_(dcd, *getMainControl()->getSelectedSystem());
			setWorkingDirFromFilename_(filename);
		}

		void DatasetControl::insertTrajectory_(TrajectoryFile* file, System& system)
			throw()
		{
			if (file->getNumberOfAtoms() != system.countAtoms())
			{
				setStatusbarText("Number of atoms do not match. Aborting...");
				delete file;
				return;
			}

			SnapShotManager* manager = new SnapShotManager(&system, 0, file);

			String name = file->getName();

			Position pos = 0;
			for (Position p = 0; p < name.size(); p++)
			{
				if (name[p] == FileSystem::PATH_SEPARATOR) pos = p;
			}
			
			if (pos) pos++;
			name = name.getSubstring(pos);
			
			QListViewItem* item = new QListViewItem(listview, name.c_str(), system.getName().c_str(), "Trajectory");
			item_to_trajectory_[item] = manager;
			insertComposite_(&system, item);
		}


		void DatasetControl::insertComposite_(Composite* composite, QListViewItem* item)
			throw()
		{
			listview->triggerUpdate();
			item_to_composite_[item] = composite;
			if (composite_to_items_.has(composite))
			{
				composite_to_items_[composite].insert(item);
			}
			else
			{
				HashSet<QListViewItem*> set;
				set.insert(item);
				composite_to_items_[composite] = set;
			}
		}

		void DatasetControl::onNotify(Message *message)
			throw()
		{
		#ifdef BALL_VIEW_DEBUG
			Log.error() << "DatasetControl "<<this<<  " onNotify " << message << std::endl;
		#endif

			GenericControl::onNotify(message);
			if (RTTI::isKindOf<RegularDataMessage>(*message))
			{
				if (RTTI::isKindOf<RegularData1DMessage>(*message))
				{
					RegularData1DMessage* ntm = RTTI::castTo<RegularData1DMessage>(*message);
					insertGrid_(ntm->getData(), (System*)ntm->getComposite(), ntm->getCompositeName());
					return;
				}
				if (RTTI::isKindOf<RegularData2DMessage>(*message))
				{
					RegularData2DMessage* ntm = RTTI::castTo<RegularData2DMessage>(*message);
					insertGrid_(ntm->getData(), (System*)ntm->getComposite(), ntm->getCompositeName());
					return;
				}
				if (RTTI::isKindOf<RegularData3DMessage>(*message))
				{
					RegularData3DMessage* ntm = RTTI::castTo<RegularData3DMessage>(*message);
					insertGrid_(ntm->getData(), (System*)ntm->getComposite(), ntm->getCompositeName());
					return;
				}
			}
			else if (RTTI::isKindOf<NewTrajectoryMessage>(*message))
			{
				NewTrajectoryMessage* ntm = RTTI::castTo<NewTrajectoryMessage>(*message);
				insertTrajectory_(ntm->getTrajectoryFile(), *(System*)ntm->getComposite());
				return;
			}
			else if (RTTI::isKindOf<CompositeMessage>(*message))
			{
				CompositeMessage *composite_message = RTTI::castTo<CompositeMessage>(*message);
				if (composite_message->getType() != CompositeMessage::REMOVED_COMPOSITE) return;
				Composite* composite = (Composite *)composite_message->getComposite();
				if (!composite_to_items_.has(composite)) return;

				// create a copy of the hashset, because it changes
				HashSet<QListViewItem*> to_delete = composite_to_items_[composite];

				HashSet<QListViewItem*>::Iterator lit = to_delete.begin();
				for (;lit != to_delete.end(); lit++)
				{
					deleteItem_(**lit);
				}
			}   
		}

		void DatasetControl::deleteItems_()
		{
			ItemList item_list = getSelectedItems(); 
			ItemList::Iterator it = item_list.begin();
			for (; it != item_list.end(); it++)
			{
				deleteItem_(**it);
			}
		}

		bool DatasetControl::deleteItem_(QListViewItem& item)
		{
			if (item_to_trajectory_.has(&item))
			{
				SnapShotManager* ssm = item_to_trajectory_[&item];
				item_to_trajectory_.erase(&item);
				delete ssm;
				setStatusbarText("deleted trajectory");
			}
			else if (item_to_grid1_.has(&item))
			{
				RegularData1D* ssm = item_to_grid1_[&item];

				RegularData1DMessage* msg = new RegularData1DMessage(RegularData1DMessage::REMOVE);
				msg->setData(*ssm);
				notify_(msg);

				item_to_grid3_.erase(&item);
				delete ssm;
				setStatusbarText("deleted 1D grid");
			}
			else if (item_to_grid2_.has(&item))
			{
				RegularData2D* ssm = item_to_grid2_[&item];

				RegularData2DMessage* msg = new RegularData2DMessage(RegularData2DMessage::REMOVE);
				msg->setData(*ssm);
				notify_(msg);

				item_to_grid2_.erase(&item);
				delete ssm;
				setStatusbarText("deleted 2D grid");
			}
			else if (item_to_grid3_.has(&item))
			{
				RegularData3D* ssm = item_to_grid3_[&item];

				RegularData3DMessage* msg = new RegularData3DMessage(RegularData3DMessage::REMOVE);
				msg->setData(*ssm);
				notify_(msg);

				item_to_grid3_.erase(&item);
				delete ssm;
				setStatusbarText("deleted 3D grid");
			}
			else
			{
				return false;
			}

			Composite* composite = item_to_composite_[&item];
			composite_to_items_[composite].erase(&item);
			item_to_composite_.erase(&item);
			GenericControl::removeItem_(&item, true);
			return true;
		}


		void DatasetControl::onContextMenu_(QListViewItem* item,  const QPoint& point, int /* column */)
		{
			if (item == 0) return;
			context_item_ = item;
			
			context_menu_.clear();
			createContextMenu_();

			insertContextMenuEntry_("Delete", SLOT(deleteItems_()));
			context_menu_.exec(point);
		}

		void DatasetControl::createContextMenu_()
		{
			if (item_to_trajectory_.has(context_item_))
			{
				insertContextMenuEntry_("Visualise/Export", SLOT(visualiseTrajectory_()));
				insertContextMenuEntry_("Load into RAM", SLOT(bufferTrajectory_()));
				return;
			}

			if (item_to_grid1_.has(context_item_) ||
		      item_to_grid2_.has(context_item_) ||
					item_to_grid3_.has(context_item_))
			{
				insertContextMenuEntry_("Save", SLOT(saveGrid_()));
				if (!item_to_grid3_.has(context_item_))
				{
					insertContextMenuEntry_("Visualise", SLOT(visualiseGrid_()));
				}
				else
				{
					insertContextMenuEntry_("ContourSurface", SLOT(computeIsoContourSurface()));
				}
			}
		}

		void DatasetControl::insertContextMenuEntry_(const QString & text, const char* member)
		{
			Index menu_entry_pos = context_menu_.insertItem(text, this, member);
			if (getSelectedItems().size() > 1) context_menu_.setItemEnabled(menu_entry_pos, false);
		}

		void DatasetControl::visualiseTrajectory_()
		{
			if (dialog_ != 0) 
			{
				dialog_->hide();
				delete dialog_;
				dialog_ = 0;
			}

			SnapShotManager* ssm = item_to_trajectory_[context_item_];

			dialog_ = new SnapshotVisualisationDialog(this);
			dialog_->setSnapShotManager(ssm);
			dialog_->show();
		}

		void DatasetControl::bufferTrajectory_()
		{
			SnapShotManager* ssm = item_to_trajectory_[context_item_];

			if (ssm->getNumberOfSnapShotsInBuffer() == 0)
			{
				if (!ssm->readFromFile())
				{
					ssm->clearBuffer();
					setStatusbarText("Could not read trajectories into buffer! Out of memory?");
				}
			}
		}

		void DatasetControl::visualiseGrid_()
		{
			if (item_to_grid1_.has(context_item_))
			{
				RegularData1D* grid = item_to_grid1_[context_item_];
				DockableRegularData1DWidget* widget = new DockableRegularData1DWidget(grid, getMainControl());
				widget->show();
				widget->zoomToFit();
				widget->undock();
			}
			else if (item_to_grid2_.has(context_item_))
			{
				RegularData2D* grid = item_to_grid2_[context_item_];
				DockableRegularData2DWidget* widget = new DockableRegularData2DWidget(grid, getMainControl());
				widget->show();
				widget->zoomToFit();
				widget->undock();
			}
		}


		void DatasetControl::saveTrajectory_()
		{
			SnapShotManager* ssm = item_to_trajectory_[context_item_];

			QString s = QFileDialog::getSaveFileName(
										getWorkingDir().c_str(),
										"DCD files(*.dcd)",
										getMainControl(),
										"Trajectory File Dialog",
										"Choose a filename to save" );

			if (s == QString::null) return;
			String filename = s.ascii();

			setWorkingDirFromFilename_(filename);

			if (!ssm->getTrajectoryFile()->copyTo(filename))
			{
				if (ssm->getTrajectoryFile()->getName() == filename)
				{
					setStatusbarText("Could not write DCDFile, you tried to save the file onto itself.", 
													 true);
				}
				else
				{
					setStatusbarText("Could not write DCDFile.", true);
				}

				return;
			}
													

			setStatusbarText("Written DCDFile", true);
		}

		String DatasetControl::chooseGridFileForOpen_()
			throw()
		{
			QString result = QFileDialog::getOpenFileName("", "*", 0, "Select a RegularData file");
			if (result == QString::null) return "";
			setWorkingDirFromFilename_(result.ascii());

			File infile;
			
			try
			{
				infile.open(result.ascii(), std::ios::in);
			}
			catch(Exception::FileNotFound)
			{
				Log.error() << "File could not be found!" << std::endl;
				return "";
			}

			return result.ascii();
		}


		void DatasetControl::add1DGrid()
			throw()
		{
			String filename = chooseGridFileForOpen_();
			if (filename == "") return;

			RegularData1D* dat = new RegularData1D;
			(*dat).binaryRead(filename);
			insertGrid_(dat, 0, filename);
			RegularData1DMessage* msg = new RegularData1DMessage(RegularData1DMessage::NEW);
			msg->setData(*dat);
			msg->setCompositeName(filename);
			notify_(msg);
		}

		void DatasetControl::add2DGrid()
			throw()
		{
			String filename = chooseGridFileForOpen_();
			if (filename == "") return;

			RegularData2D* dat = new RegularData2D;
			(*dat).binaryRead(filename);
			insertGrid_(dat, 0, filename);
			RegularData2DMessage* msg = new RegularData2DMessage(RegularData2DMessage::NEW);
			msg->setData(*dat);
			msg->setCompositeName(filename);
			notify_(msg);
		}


		void DatasetControl::add3DGrid()
			throw()
		{
			String filename = chooseGridFileForOpen_();
			if (filename == "") return;

			RegularData3D* dat = new RegularData3D;
			(*dat).binaryRead(filename);
			insertGrid_(dat, 0, filename);
			RegularData3DMessage* msg = new RegularData3DMessage(RegularData3DMessage::NEW);
			msg->setData(*dat);
			msg->setCompositeName(filename);
			notify_(msg);
		}

		void DatasetControl::insertGrid_(RegularData1D* data, System* system, const String& name)
			throw()
		{
			QListViewItem* item = createListViewItem_(system, name, "1D Grid");
			item_to_grid1_[item] = data;
		}

		void DatasetControl::insertGrid_(RegularData2D* data, System* system, const String& name)
			throw()
		{
			QListViewItem* item = createListViewItem_(system, name, "2D Grid");
			item_to_grid2_[item] = data;
		}

		void DatasetControl::insertGrid_(RegularData3D* data, System* system, const String& name)
			throw()
		{
			QListViewItem* item = createListViewItem_(system, name, "3D Grid");
			item_to_grid3_[item] = data;
		}


		QListViewItem* DatasetControl::createListViewItem_(System* system, const String& name, const String& type)
			throw()
		{
			QListViewItem* item;
			if (system != 0) 
			{
				item = new QListViewItem(listview, name.c_str(), system->getName().c_str(), type.c_str());
				insertComposite_(system, item);
			}
			else
			{ 	
				item = new QListViewItem(listview, name.c_str(), "", type.c_str());
			}

			return item;
		}


		String DatasetControl::chooseGridFileForSave_()
			throw()
		{
			QString qs = QFileDialog::getSaveFileName("", "*", 0, "Select a RegularData file");
			if (qs == QString::null) return "";
			setWorkingDirFromFilename_(qs.ascii());

			String result = qs.ascii();
			if (result.isEmpty()) return 0;

			File test;
			
			try
			{
				test.open(result, std::ios::out);
			}
			catch(Exception::GeneralException)
			{
				Log.error() << "File could not be written!" << std::endl;
				return "";
			}

			return result;
		}

		void DatasetControl::saveGrid_()
			throw()
		{ 
			String filename = chooseGridFileForSave_();
			if (filename == "") return;

			if (item_to_grid1_.has(context_item_))
			{
				item_to_grid1_[context_item_]->binaryWrite(filename);
			}
			else if (item_to_grid2_.has(context_item_))
			{
				item_to_grid2_[context_item_]->binaryWrite(filename);
			}
			else if (item_to_grid3_.has(context_item_))
			{
				item_to_grid3_[context_item_]->binaryWrite(filename);
			}

			setStatusbarText("Grid successfully written...");
		}

		void DatasetControl::updateSelection()
			throw()
		{
			GenericControl::updateSelection();

			QListViewItemIterator it(listview);
			for (; it.current(); ++it)
			{
				QListViewItem* item = it.current();
				if (!item->isSelected()) continue;

				if (item_to_grid1_.has(item))
				{
					RegularData1DMessage* message = new RegularData1DMessage(RegularData1DMessage::SELECTED);
					message->setData(*item_to_grid1_[item]);
					message->setCompositeName(item->text(0).ascii());
					notify_(message);
					break;
				}

				if (item_to_grid2_.has(item))
				{
					RegularData2DMessage* message = new RegularData2DMessage(RegularData2DMessage::SELECTED);
					message->setData(*item_to_grid2_[item]);
					message->setCompositeName(item->text(0).ascii());
					notify_(message);
					break;
				}

				if (item_to_grid3_.has(item))
				{
					RegularData3DMessage* message = new RegularData3DMessage(RegularData3DMessage::SELECTED);
					message->setData(*item_to_grid3_[item]);
					message->setCompositeName(item->text(0).ascii());
					notify_(message);
					break;
				}
			}
		}

		List<std::pair<RegularData3D*, String> > DatasetControl::get3DGrids()
			throw()
		{
			List<std::pair<RegularData3D*,String> > grids;
			HashMap<QListViewItem*, RegularData3D*>::Iterator it = item_to_grid3_.begin();
			for (; it != item_to_grid3_.end(); it++)
			{
				std::pair<RegularData3D*, String> p((*it).second, (*it).first->text(0).ascii());
				grids.push_back(p);
			}
			return grids;
		}


		void DatasetControl::computeIsoContourSurface()
		{
			// execute the surface dialog and abort if cancel is clicked
			if (surface_dialog_ == 0) 
			{
				surface_dialog_ = new ContourSurfaceDialog(this, "ContourSurfaceDialog");
				surface_dialog_->setDatasetControl(this);
				registerWidgetForHelpSystem(surface_dialog_, "datasetControl.html#isocontour_surfaces");
			}
			if (!surface_dialog_->exec()) return;

			// Create a new contour surface.
			ContourSurface cs(*surface_dialog_->getGrid(), surface_dialog_->getThreshold());

			if (cs.vertex.size() == 0)
			{
				setStatusbarText("Could not calculate ContourSurface, no grid points found for threshold!", true);
				surface_dialog_->show();
				return;
			}

			Mesh* mesh = new Mesh;
			
			// reorient the vertices for OpenGL
			for (Position t = 0; t < cs.triangle.size(); t++)
			{
				const Mesh::Triangle& org = cs.triangle[t];
				Mesh::Triangle tri;
				tri.v1 = org.v1;
				tri.v2 = org.v3;
				tri.v3 = org.v2;
				mesh->triangle.push_back(tri);
			}

			mesh->vertex = cs.vertex;
			mesh->normal = cs.normal;

			mesh->colors.clear(); 
			mesh->colors.push_back(surface_dialog_->getColor());

			//////////////////////////////////////////////
			// Create a new representation containing the contour surface.
			Representation* rep = new Representation();
			rep->insert(*mesh);
			rep->setModelType(MODEL_CONTOUR_SURFACE); 

			// Make sure BALLView knows about the new representation.
			getMainControl()->insert(*rep);
			getMainControl()->update(*rep);
		}


		DatasetControl::DatasetControl(const DatasetControl& control)
			throw()
			: GenericControl(control)
		{
		}

	} // namespace VIEW
} // namespace BALL
