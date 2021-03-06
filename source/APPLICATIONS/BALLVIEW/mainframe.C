// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// $Id: mainframe.C,v 1.55.2.19 2005/10/25 14:16:13 amoll Exp $
//

#include "mainframe.h"
#include "aboutDialog.h"
#include "icons.h"
#include "demoTutorialDialog.h"

#include <BALL/VIEW/KERNEL/moleculeObjectCreator.h>
#include <BALL/VIEW/KERNEL/server.h>
#include <BALL/VIEW/RENDERING/POVRenderer.h>
#include <BALL/VIEW/RENDERING/VRMLRenderer.h>
#include <BALL/VIEW/WIDGETS/molecularStructure.h>
#include <BALL/VIEW/WIDGETS/molecularControl.h>
#include <BALL/VIEW/WIDGETS/geometricControl.h>
#include <BALL/VIEW/WIDGETS/logView.h>
#include <BALL/VIEW/WIDGETS/helpViewer.h>
#include <BALL/VIEW/DIALOGS/downloadPDBFile.h>
#include <BALL/VIEW/DIALOGS/labelDialog.h>
#include <BALL/VIEW/DIALOGS/displayProperties.h>

#include <BALL/FORMAT/PDBFile.h>
#include <BALL/COMMON/version.h>

#ifdef BALL_PYTHON_SUPPORT
#	include <BALL/VIEW/WIDGETS/pyWidget.h>
#endif

#include <qmenubar.h>
#include <qlabel.h>
#include <qprinter.h>
#include <qpainter.h>
#include <qimage.h>
#include <qmessagebox.h>
#include <qcursor.h>

#include <sstream>

namespace BALL
{
	using namespace std;
	using namespace BALL::VIEW;

	Mainframe::Mainframe(QWidget* parent, const char* name)
		:	MainControl(parent, name, ".BALLView"),
			scene_(0),
			fullscreen_(false)
	{
		#ifdef BALL_VIEW_DEBUG
		Log.error() << "new Mainframe " << this << std::endl;
		#endif

		// ---------------------
		// setup main window
		// ---------------------
		setCaption("BALLView");
		setIcon(QPixmap(bucky_64x64_xpm));
		resize(800,600);
		// make sure submenus are the first 
		initPopupMenu(FILE_OPEN);
		initPopupMenu(EDIT);
		initPopupMenu(BUILD);
		initPopupMenu(DISPLAY);
		initPopupMenu(DISPLAY_VIEWPOINT);
		initPopupMenu(MOLECULARMECHANICS);
		initPopupMenu(TOOLS);
	#ifdef BALL_PYTHON_SUPPORT
		initPopupMenu(TOOLS_PYTHON);
	#endif
		initPopupMenu(WINDOWS);

		// ---------------------
		// Logstream setup -----
		// ---------------------
//		Lo g.remove(std::cout);
// 		Log.remove(std::cerr);
		setLoggingFilename("BALLView.log");
		setAcceptDrops(true);

		CHECK_PTR(new MolecularControl(		this, "Structures"));
		CHECK_PTR(new GeometricControl(		this, "Representations"));
		CHECK_PTR(new DatasetControl(			this, "Datasets"));
		CHECK_PTR(new DisplayProperties(	this, "DisplayProperties"));
		CHECK_PTR(new LabelDialog(				this, "LabelDialog"));
		CHECK_PTR(new MolecularFileDialog(this, "MolecularFileDialog"));
		CHECK_PTR(new DownloadPDBFile(		this, "DownloadPDBFile", false));
		CHECK_PTR(new MolecularStructure(	this, "MolecularStructure"));
		CHECK_PTR(new HelpViewer(					this, "Documentation"));
		CHECK_PTR(new LogView(						this, "Logs"));

		Scene::stereoBufferSupportTest();
		scene_ = new Scene(this, "3D View");
		CHECK_PTR(scene_);
		setCentralWidget(scene_);

		// setup the VIEW server
		Server* server = new Server(this);
		CHECK_PTR(server);
		// registering object generator
		MoleculeObjectCreator* object_creator = new MoleculeObjectCreator;
		server->registerObjectCreator(*object_creator);

		DemoTutorialDialog* demo = new DemoTutorialDialog(this, "BALLViewDemo");
		CHECK_PTR(demo);

		#ifdef BALL_PYTHON_SUPPORT
			new PyWidget(this, "Python Interpreter");
		#endif

		// ---------------------
		// Menus ---------------
		// ---------------------
		String hint;

		// File Menu
		Index entry = insertMenuEntry(MainControl::FILE_EXPORT, "POVRa&y scene", this, SLOT(exportPOVRay()), CTRL+Key_Y);

		// registerMenuEntryForHelpSystem
		RegisterHelpSystemMessage* msg = new RegisterHelpSystemMessage();
		msg->setMenuEntry(entry);
		msg->setURL("tips.html#povray");
		notify_(msg);

		insertMenuEntry(MainControl::FILE, "Print", this, SLOT(printScene()));

		insertMenuEntry(MainControl::FILE_OPEN, "Project", this, 
										SLOT(loadBALLViewProjectFile()), 0, 2);

		save_project_id_ = insertMenuEntry(MainControl::FILE, "Save Project", this, 
										SLOT(saveBALLViewProjectFile()));

		// Display Menu
		insertMenuEntry(MainControl::DISPLAY, "Toggle Fullscreen", this, SLOT(toggleFullScreen()),
										ALT+Key_X);

		// Help-Menu -------------------------------------------------------------------
		insertMenuEntry(MainControl::HELP, "Demo", demo, SLOT(showDemo()));
		insertMenuEntry(MainControl::HELP, "Tutorial", demo, SLOT(showTutorial()));
		insertPopupMenuSeparator(MainControl::HELP);
		insertMenuEntry(MainControl::HELP, "About", this, SLOT(about()));

		// Menu ------------------------------------------------------------------------
		menuBar()->setSeparator(QMenuBar::InWindowsStyle);

//   		#ifdef BALL_QT_HAS_THREADS
			stop_simulation_id_ = insertMenuEntry(MainControl::MOLECULARMECHANICS, "Abort Calculation", this, 
											SLOT(stopSimulation()), ALT+Key_C);
			insertPopupMenuSeparator(MainControl::MOLECULARMECHANICS);
			setMenuHint(stop_simulation_id_, "Abort a running simulation thread");
//   		#endif

		complement_selection_id_ = insertMenuEntry(MainControl::EDIT, "Toggle Selection", this, SLOT(complementSelection()));

		setStatusbarText("Ready.");
	}

	Mainframe::~Mainframe()
		throw()
	{
	}

	void Mainframe::exportPOVRay()
	{
		QString qresult = QFileDialog::getSaveFileName(
											getWorkingDir().c_str(),
											"*.pov",
											this,
											"Export POVRay File",
											"Choose a filename" );

	 	if (qresult == QString::null) return;

		String result = qresult.ascii();

		if (!result.hasSuffix(".pov"))
		{
			result += ".pov";
		}

		POVRenderer pr(result);
		if (!scene_->exportScene(pr))
		{
			setStatusbarText("Could not export POV to " + result, true);
		}
		else
		{
			setStatusbarText("Exported POV to " + result);
		}
	}

	void Mainframe::toggleFullScreen()
	{
		if (!fullscreen_)
		{
			// This call is needed because showFullScreen won't work
			// correctly if the widget already considers itself to be fullscreen.
			showNormal();	
			showFullScreen();
		}
		else
		{
			showNormal();
		}
		fullscreen_ = !fullscreen_;
	}

	void Mainframe::openFile(const String& file)
		throw()
	{
		if (composites_locked_) return;

		setStatusbarText(String("Opening file ") + file + "...");

		if (file.hasSuffix(".bvp"))
		{
			MainControl::loadBALLViewProjectFile(file);
			return;
		}

#ifdef BALL_PYTHON_SUPPORT
		if (file.hasSuffix(".py"))
		{
			PyWidget::getInstance(0)->run(file);
			return;
		}
#endif

		if (file.hasSuffix(".dcd"))
		{
			DatasetControl::getInstance(0)->addTrajectory(file);
			return;
		}

		MolecularFileDialog::getInstance(0)->openFile(file);
	}


	void Mainframe::printScene()
	{
		QPrinter printer;
		if (!printer.setup(this)) return;
		
		setStatusbarText("printing..");

		QPainter p;
		if(!p.begin(&printer)) return; 

		QImage pic = scene_->grabFrameBuffer();
		p.drawImage(0,0, pic);	
		p.end();

		setStatusbarText("finished printing");
	}

	
	void Mainframe::saveBALLViewProjectFile()
		throw()
	{
		QString qresult = QFileDialog::getSaveFileName(
											getWorkingDir().c_str(), 
											"*.bvp", 
											this, 
											"Select a BALLView project file");

	 	if (qresult == QString::null) return;

		String result = qresult.ascii();
		if (result.isEmpty()) return;

		if (!result.hasSuffix(".bvp")) result += ".bvp";

 		MainControl::saveBALLViewProjectFile(result);
		setStatusbarText("Saved project to " + result);
	} 


	void Mainframe::loadBALLViewProjectFile()
		throw()
	{
		QString result = QFileDialog::getOpenFileName(
				getWorkingDir().c_str(), "*.bvp", 0, "Select a BALLView project file");
		if (result.isEmpty())
		{
			return;
		}

 		MainControl::loadBALLViewProjectFile(result.ascii());
	}


	void Mainframe::keyPressEvent(QKeyEvent* e)
	{
		if (e->key() == Key_Escape) 
		{
			Scene::getInstance(0)->switchToLastMode();
			return;
		}

		if (e->key() == Key_Enter) 
		{
			if (composite_manager_.getNumberOfComposites() == 0) return;

			if (getMolecularControlSelection().size() == 0)
			{
				control_selection_.push_back(*composite_manager_.begin());
			}
				
			MolecularStructure::getInstance(0)->centerCamera();
			return;
		}

		// check all menu entries if Alt or CTRL is pressed to enable shortcuts
		if (e->key() == Qt::Key_Alt ||
				e->key() == Qt::Key_Control)				
		{
			checkMenus();
			return;
		}

		#ifdef BALL_PYTHON_SUPPORT
			PyWidget::getInstance(0)->reactTo(*e);
			e->accept();
		#endif
	}

	
	void Mainframe::reset()
	{
		if (composites_locked_ || getPrimitiveManager().updateRunning()) return;

		DisplayProperties* dp = DisplayProperties::getInstance(0);
		dp->setDrawingPrecision(DRAWING_PRECISION_HIGH);
		dp->selectModel(MODEL_STICK);
		dp->selectColoringMethod(COLORING_ELEMENT);
		dp->selectMode(DRAWING_MODE_SOLID);
		dp->setTransparency(0);
		dp->setSurfaceDrawingPrecision(6.5);

		// remove all loaded Composites
		HashSet<Composite*> composites = getCompositeManager().getComposites();
		HashSet<Composite*>::Iterator it = composites.begin();
		
		for (; +it; ++it)
		{
			remove(**it, true, false);
		}

		// remove all Representations
		PrimitiveManager::RepresentationList reps = getPrimitiveManager().getRepresentations();
		PrimitiveManager::RepresentationList::Iterator rit = reps.begin();

		for (; rit != reps.end(); ++rit)
		{
			remove(**rit);
		}
	}

	
	void Mainframe::checkMenus()
	{
		MainControl::checkMenus();
		menuBar()->setItemEnabled(save_project_id_, !composites_locked_);
	}

	void Mainframe::about()
	{
		// Display about dialog
		AboutDialog about;
		QString version = QString("QT ") + qVersion();
#ifdef BALL_QT_HAS_THREADS
		version += "(mt)";
#endif
		about.qt_version_label->setText(version);
		QFont font = about.BALLView_version_label->font();
		about.BALLView_version_label->setText(QString("BALLView ") + BALL_RELEASE_STRING);
		font.setPixelSize(18);
		about.BALLView_version_label->setFont(font);
		about.BALL_version_label->setText(__DATE__);
		about.exec(); 
	}

}
