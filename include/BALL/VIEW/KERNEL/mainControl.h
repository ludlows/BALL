// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// $Id: mainControl.h,v 1.72.2.23 2005/11/26 15:49:13 amoll Exp $
//

#ifndef BALL_VIEW_KERNEL_MAINCONTROL_H
#define BALL_VIEW_KERNEL_MAINCONTROL_H

#ifndef BALL_CONCEPT_EMBEDDABLE_H
#	include <BALL/CONCEPT/embeddable.h>
#endif

#ifndef BALL_DATATYPE_HASHMAP_H
#	include <BALL/DATATYPE/hashMap.h>
#endif

#ifndef BALL_VIEW_KERNEL_CONNECTIONOBJECT_H
#	include <BALL/VIEW/KERNEL/connectionObject.h>
#endif

#ifndef BALL_VIEW_KERNEL_PRIMITIVEMANAGER_H
#	include <BALL/VIEW/KERNEL/primitiveManager.h>
#endif

#ifndef BALL_VIEW_KERNEL_COMPOSITEMANAGER_H
#	include <BALL/VIEW/KERNEL/compositeManager.h>
#endif

#ifndef BALL_FORMAT_INIFILE_H
# include <BALL/FORMAT/INIFile.h>
#endif

#ifndef BALL_SYSTEM_FILE_H
# include <BALL/SYSTEM/file.h>
#endif

#ifndef BALL_STRUCTURE_FRAGMENTDB_H
# include <BALL/STRUCTURE/fragmentDB.h>
#endif

#include <qmainwindow.h>
#include <qapplication.h>
#include <qmenubar.h>    // menus
#include <qlabel.h>			 // statusbar
#include <qtimer.h>

namespace BALL
{
	namespace VIEW
	{
		class ModularWidget;
		class Preferences;
		class MainControlPreferences;
		class NetworkPreferences;
		class GeometricObjectSelectionMessage;
		class SimulationThread;

		/**	MainControl is the main administration unit for a program and must be
				used by all	applications.
				It is a storage facility for Composite objects, the graphical
				Representation and the inserted ModularWidget.
				The interface for the Composite administration is implemented in CompositeManager, and for
				the Representation 's in PrimitiveManager.
				\par
				This class is also the root ConnectionObject. 
				Therefore all messages will be handled from this class. 
				\par
				It is also derived from the QT::QMainWindow and therefore the main
				widget of an application must be derived from this class. Further it has the
				necessary interface methods to create and update the menus of the main application.
				\par
				It handles also the general preferences tab Preferences of the main application and notifies all
				registered ModularWidget objects if the preferences have changed. \par
				The preferences of the application are stored in an INIFile.
				The default name of this file is ".BALLView".\par
				\par
				<b>Caveat:</b> Due to a peculiarity of the QT Meta Object Compiler (MOC)
				you have to specify the full namespace qualified name of this class when deriving from it. \par
				So don't use\par 
				<tt> class foo : public MainControl </tt>; but \par
				<tt> class foo : public BALL::VIEW::MainControl </tt> instead. 
		\ingroup ViewKernelConnectivity
		*/
		class BALL_VIEW_EXPORT MainControl
			: public QMainWindow,
				public ConnectionObject,
				public Embeddable
		{
			friend class PrimitiveManager;
			friend class SimulationThread;

			Q_OBJECT

			public:

			BALL_EMBEDDABLE(MainControl,Embeddable)

			/**	@name Enumerations
			*/
			//@{
			/**	Standard Popup menu IDs.
					This enum defines symbolic names for the menu IDs of the most common popup menus.
					The popups are created, if a menu entry is requested for any of the popups.
					\see insertMenuEntry
			*/
			enum PopUpID
			{
				/// File menu
				FILE        = 10001,

				/// File menu sub menu open
				FILE_OPEN,

				/// File menu sub menu import [currently unused]
				FILE_IMPORT,

				/// File menu sub menu export
				FILE_EXPORT,

				/// Edit menu
				EDIT = 10100,

				/// Build menu 
				BUILD = 10200,

				/// Display menu 
				DISPLAY = 10300,

				/// Display Viewpoint submenu
				DISPLAY_VIEWPOINT,
				
				/// Display Stereo submenu
				DISPLAY_STEREO,
				
				/// Display Animation submenu
				DISPLAY_ANIMATION,

				/// Simulations menu
				MOLECULARMECHANICS = 10400,
				
				/// Molmec submenu for force field selection
				CHOOSE_FF,

				/// Tools menu
				TOOLS = 10500,

				/// Create grid submenu in Tools
				TOOLS_CREATE_GRID,

				/// Python submenu in Tools
				TOOLS_PYTHON,

				/// Windows menu
				WINDOWS = 10600,

				/// Userdefined menus
				USER = 10700,

				/// Help menu
				HELP = 10800
			};

			//@}
			/**	@name	Constructors and Destructor
			*/
			//@{

			/** Default Constructor.
					Reads the the INIFile <tt>inifile</tt> and connects the qt signal 
					<b> aboutToQuit </b> with the slot aboutToExit().
					The state of the MainControl is:
						-  no Composite objects stored
						-  no Representation objects stored
						-  no general Preferences dialog added
						-  no MainControlPreferences dialog added
					\par
					\param  parent the new parent widget
					\param  name the new name of the widget
					\param  inifile the new preferences filename
					\see    INIFile
					\see    Preferences
			*/
			MainControl(QWidget* parent = 0, const char* name = 0 , String inifile = ".BALL.preferences")
				throw();

			/** Destructor.
					Calls clear
			*/
			virtual ~MainControl()
				throw();
					
			// copy ctor needed for Python support only!
			MainControl(const MainControl& main_control)
				throw();

			/** Explicit default initialization.
			*/
			virtual void clear()
				throw();

			//@}
			/**	@name	Accessors
			*/
			//@{

			/** Get the primitive manager.
					The class PrimitiveManager contains all Representation objects and GeometricObject.
			*/
			PrimitiveManager& getPrimitiveManager()
				throw() { return primitive_manager_;}

			/** Get the composite manager.
					The class CompositeManager is the owner of all Composite objects.
			*/
			CompositeManager& getCompositeManager()
				throw() { return composite_manager_;}
			
			/** Redraws all Representation objects for a Composite.
					If the Composite is not inserted into this MainControl <tt>false</tt> will be returned.
					updateRepresentationsOf() is called after receiving a CompositeMessage with type CHANGED_COMPOSITE in onNotify().
					It sends a RepresentationMessage with type UPDATE for every Representation, which was build for the 
					Composite.  After this a SceneMessage is send to redraw the Scene.
					Remember:
					If you changed a Composite in MainControl or a derived class, the MainControl doesnt get 
					notified, from the CompositeMessage it sends. So you have to call this function instead 
					of sending the message.
					\param  composite the Composite that should be updated
					\param  rebuild if set to true, the model is rebuilded, otherwise just the coloring is updated
					\param  force is set to true, also rebuild non surface models (only usefull with rebuild = true)
					\return true if an update was performed
			*/
			bool updateRepresentationsOf(const Composite& composite, bool rebuild = true, bool force = false)
				throw();

			/** Redraws all inserted Representation, but doesnt change the Models.
					\param rebuild_display_lists set to true lets the Scene rebuild the GLDisplayList objects.
					\see updateRepresentationsOf
			*/
			void redrawAllRepresentations(bool rebuild_display_lists = false)
				throw();

			/** Update a Composite in all ModularWidget.
			 		This method differs wheter the composites hierarchy was changed or not.
					The update is faster if the hierarchy is unchanged, because e.g. the
					MolecularControl doesnt have to rebuild the ListViewItem tree.
			 		A CompositeMessage with type CHANGED_COMPOSITE or CHANGED_COMPOSITE_HIERARCHY is send and
					updateRepresentationsOf(composite) is called.
					\return false if the CompositeManager doesnt contain the Composite
			*/
			void update(Composite& composite, bool changed_hierarchy = true)
				throw();

			/** Insert a Composite and notify all ModularWidget.
			 		The Composite has to be created on the heap!!!
			 		A CompositeMessage with type NEW_COMPOSITE is send and
					CompositeManager::insert called.
					\return false if the CompositeManager contains the Composite
			*/
			bool insert(Composite& composite, String name = "")
				throw();

			/** Remove a Composite and notify all ModularWidget.
			 		A CompositeMessage with type REMOVED_COMPOSITE is send and
					CompositeManager::remove called.
					@param update update Representations if needed
					\return false if the CompositeManager doesnt contain the Composite
			*/
			bool remove(Composite& composite, bool to_delete = true, bool update = true)
				throw();

			/** Update a Representation
			 		A RepresentationMessage with type UPDATE and a SceneMessage is send.
					\return false if the PrimitiveManager doesnt contain the Representation
			*/
			bool update(Representation& rep)
				throw();

			/** Insert a Representation
			 		The Representation must be created on the heap!!!
			 		A RepresentationMessage with type NEW is send.
					\return false if the PrimitiveManager contains the Representation
			*/
			bool insert(Representation& rep)
				throw();

			/** Remove a Representation
			 		A RepresentationMessage with type REMOVE is send.
					\return false if the PrimitiveManager doesnt contain the Representation
			*/
			bool remove(Representation& rep)
				throw();

				
			/** Mutable inspection of the preferences dialog.
					\return   Preferences* a pointer to the Preferences dialog, (<tt> 0</tt> if not present)
			*/
			Preferences* getPreferences()
				throw();

			/** Mutable inspection of the INIFile.
			*/
			INIFile& getINIFile()
				throw();

			/** Non-mutable inspection of the INIFile.
			*/
			const INIFile& getINIFile() const
				throw();

			/** Message handling method.
					Handles messages sent by other registered ModularWidget objects.
					Virtual function for overriden the message handling system. 
					Take care to call this function in your own ModularWidget.
					There is no need to call this function, because it will be called from the 
					message handling mechanism.
					<b>Remember:</b> A ModularWidget is not notified by the Messages it sends itself!
					\param message the pointer to the message that should be processed
					\see   ModularWidget
					\see   Message
			*/
			virtual void onNotify(Message *message)
				throw();

			/** Send a Message from Python.
					Otherwise, you should prefer to use ModularWidget::notify_.
					The MainControl itself also reacts to a Message, send with this method.
					The Message will be deleted, after it was send to all ModularWidget's.
			*/
			void sendMessage(Message& message)
				throw();

			/// Get the ID of the last highlighted menu entry (used for the HelpViewer)
			Index getLastHighLightedMenuEntry() { return last_highlighted_menu_entry_;}
			
			public slots:

			//@}
			/** @name Public slots
			*/
			//@{

			/**	Initialize all registered ModularWidget objects.
					It initializes the menu structure, the preferences dialogs and connects
					every ModularWidget with the MainControl.
					This method also creates the first menu entry <b> FILE </b> with its subentry
					<b> EXIT </b> to exit the application.
					See ModularWidget for further information concerning menu structure creation
					and preferences handling. \par
					Calls registerConnectionObject() \par
					Calls fetchPreferences() \par
					Calls applyPreferences() \par
					Calls insertMenuEntry() \par
					Calls ModularWidget::initializeWidget() \par
					Calls QMainWindow::show() \par
					Note: Call this method to start the application.
			*/
			virtual void show();

			/** Menu checking method.
					This method checks, enables or disables all inserted menu entries of the
					registered ModularWidget objects.
					If this method is overridden make sure that it will be called at the end of 
					the new <b>checkMenus()</b> method.
					See ModularWidget for further information concerning menu structure creation.\par
					<b>Note:</b> This method will be called internally whenever the menu structure needs an update.
					Calls ModularWidget::checkMenu\par
					\see        ModularWidget::checkMenu
			*/
			virtual void checkMenus();

			/// Stop a currently running calculation
			void stopSimulation();
	
			///
			void complementSelection();

			/** Last second cleanup.
					This method will be called internally if the MainControl is about to be destroyed.
					This method stores the preferences and finalizes all ModularWidget objects
					and the MainControl.
					Must be called after your own cleanup routine if you override this method.\par
					Calls ModularWidget::finalizePreferencesTab \par
					Calls ModularWidget::finalizeWidget \par
					Calls writePreferences \par
					Calls finalizePreferencesTab \par
					Calls removeModularWidget \par
					Calls INIFile::write \par
			*/
			virtual void aboutToExit();

			/** Slot that is called when a menu item is highlighted.
					It is used to show a hint for every menu entry.
					@see setMenuHint
					@see getMenuHint
			*/
			void menuItemHighlighted(int id)
				throw();
			
			/// Interface to QT events, e.g. to communicate with other threads
			virtual void customEvent( QCustomEvent * e );

			/// Make the program exit
			virtual void quit();

			///
			bool isAboutToQuit() { return about_to_quit_;}
			
			/// overloaded from QT for Python Interface
			virtual void resize (int w, int h );

			public:
			
			//@}
			/**	@name	Automatic module registration, menu construction and preferences handling.
			*/
			//@{

			/** Return the MainControl of an QObject.
					This method returns the MainControl that should be the root of the
					ConnectionObject tree from a given widget or dialog QObject.
					Because we use the qt library, every widget or dialog has 
					QObject as a base class. MainControl is the main application and therefore
					all widgets and dialogs are its children. We use the qt
					QObject tree mechanism to return the MainControl for a given QObject.\par
					<b>Note</b>: This method is used internally from the ModularWidget registration process.
					\return   MainControl* the root of the ConnectionObject tree
					\see      ConnectionObject
					\see      ModularWidget
			*/
			static MainControl* getMainControl(const QObject* object)
				throw();
			
			/** Insert a new menu entry into menu <b>ID</b> 
					(creates a new menu if <b>ID</b> not existent).
					See the documentation of the qt library for more information concerning menu creation.
					\param ID the menu ID to which the new menu entry should be inserted
					\param name the name of the new menu entry
					\param receiver the object to which the menu action will be connected
					\param slot the function that will be called by activation of the menu entry
					\param accel the acceleration key
					\param entry_ID the id for the new menu entry (default: -1, will create a new one)
					\param hint
					\return int the new entry_ID
			*/
			Index insertMenuEntry (Index parent_id, const String& name, const QObject* receiver = 0, 
													 const char* slot = 0, Index accel = 0, Index pos = -1)
				throw();

			/// 
			void removeMenuEntry (Index parent_id, Index entry_ID)
				throw();
			
			/**	Initialize a new popup menu <b> ID</b>. 
					If the MainControl has already the popup menu <b>ID</b> that QPopupMenu is returned.
					See the documentation of the qt library for more information concerning the class QPopupMenu.
					\param    ID the ID of the menu entry to be created.
					\return   QPopupMenu* a pointer to the created QPopupMenu
					\see      PopUpID
			*/	
			virtual QPopupMenu* initPopupMenu(int ID)
				throw();

			/** Insert a separator into the popup menu <b> ID</b>. 
					If the menu <b>ID</b> is not existent, it will be created first.
					\param ID the id of the menu to which a separator will be inserted
					\see   PopUpID
			*/
			void insertPopupMenuSeparator(int ID)
				throw();
			
			/** Apply all preferences.
					This method is called automatically by applyPreferencesClicked() and calls
					applyPreferences() for all registered ModularWidgets.
					<b>Note:</b> If this method is overridden, call this method at the end of the
					overriden method to make sure that the general preferences are applied.
					\see    ModularWidget
					\see    Preferences
			*/
			virtual void applyPreferences()
				throw();
			
			/** Fetch the preferences from the INIfile.
					Calls fetchPreferences() for all registered ModularWidgets.
					<b>Note:</b>If this method is overridden, call this method at the end of the
					overriden method to make sure that the general preferences are fetched.
					\param  inifile the INIFile that contains the needed values
			*/
			virtual void fetchPreferences(INIFile &inifile)
				throw();
			
			/** Writes the widgets preferences to the INIFile.
					Calls writePreferences() for all registered ModularWidgets and
					Preferences::savePreferences().
					<b>Note:</b> If this method is overridden, call this method at the end of the
					overriden method to make sure that the general preferences are written.
					\param  inifile the INIFile that contains the needed values
			*/
			virtual void writePreferences(INIFile &inifile)
				throw();
			
			/// Restore the positions the main window and of all DockWindow's from the INIFile
			virtual void restoreWindows(const INIFile& inifile)
				throw();
			
			/** Add a new ModularWidget to this MainControl.
					This method will be called internally by the ModularWidget registration process.
					\param  widget the ModularWidget to be inserted into this mainControl
			*/
			void addModularWidget(ModularWidget* widget)
				throw();

			/** Remove a ModularWidget from the MainControl.
					This method will be called internally by the ModularWidget registration process.
					\param  widget the ModularWidget to be removed
			*/
			void removeModularWidget(ModularWidget* widget)
				throw();

			//@}
			/**	@name	Accessors and Settings
			*/
			//@{
			
			/// Get the HashSet with the selected (e.g. picked) Composite objects (const)
			const HashSet<Composite*>& getSelection() const
				throw();

			/// Get the HashSet with the selected (e.g. picked) Composite objects
			HashSet<Composite*>& getSelection() 
				throw();

			/// Get the selection (highlighted items) of the MolecularControl (not the selection with checkboxes)
			List<Composite*>& getMolecularControlSelection()
				throw();

			/// If exactly one System is selected in the Control, return a pointer to this system, otherwise 0.
			System* getSelectedSystem()
				throw();

			///	Select a Composite recursive and add all Atom and AtomContainer objects to the selection.
			void selectCompositeRecursive(Composite* composite, bool first_call=false)
				throw();

			/// Select a Composite recursive and add all Atom and AtomContainer objects to the selection.
			void deselectCompositeRecursive(Composite* composite, bool first_call=false)
				throw();

			/** Clear Selection
			 		Deselect all Composites and clear the selection list in the MainControl
			*/
			void clearSelection()
				throw();

			/** Print some informations for the selection in the statusbar.
					Called by selectComposites_().
					If one Atom is selected, its position is printed.
					If two Atom objects are selected, their distance,
					for three Atom 's their angle and
					for four Atom 's their torsion angle.
					Else the number of items is printed.
			*/
			void printSelectionInfos()
				throw();

			/** Sets the text in the statusbar.
					The statusbar has a label, whose text is set to the given argument.
					It is possible to notify the user with a beep sound.
			*/
			void setStatusbarText(const String& text, bool important = false, bool beep = false)
				throw();

			///
			String getStatusbarText() const
				throw();

			/// Set a hint for a menu entry
			void setMenuHint(Index id, const String& hint)
				throw();

			/// Get the hint for a menu entry
			const String& getMenuHint(Index id) const
				throw();
			
			/// Get a const reference for the fragment database
			const FragmentDB& getFragmentDB() const
				throw() { return fragment_db_;}

			/** Check wheter the stored composites can be modified at the moment.
					This method returns true e.g. while a MD simulation is running.
			*/
			bool compositesAreLocked() throw();

			///
			bool lockCompositesFor(ModularWidget* widget) throw();

			///
			bool unlockCompositesFor(ModularWidget* widget) throw();

			///
			ModularWidget* getLockingWidget() throw();

			///
			bool updateOfRepresentationRunning() throw();
					
			/// Returns true, if the simulation was told to stop, but hasnt done this so far.
			bool stopedSimulation() { return stop_simulation_;}

			/** Set the simulation thread.
			 		The instance of SimulationThread will be deleted after it
					has finished. If an other simulation is still running, this
					method returns false.
			*/
			bool setSimulationThread(SimulationThread* thread)
				throw();

			/** Get the currently running SimulationThread or
			 		zero pointer if no simulation running.
			*/
			SimulationThread* getSimulationThread()
				throw();

			/** Enable the delete entry for GenericControls.
					Called by a GenericControl, if it has a selection, that can be deleted.
			*/
			void setDeleteEntryEnabled(bool state)
				throw();
	
			/** Insert the delete entry for GenericControls.
					Called by all GenericControls.
			*/
			void insertDeleteEntry()
				throw();

			///
			String getWorkingDir() const
				throw() { return working_dir_;}

			///
			void setWorkingDir(const String& dir)
				throw();

			///
			void enableLoggingToFile()
				throw();
			
			///
			void disableLoggingToFile()
				throw();

			///
			void setLoggingFilename(const String& string)
				throw();

			///
			const String& getLoggingFilename() const
				throw();

			///
			void setProxy(const String& host, Position port);

			///
			String getProxy() const { return proxy_;}

			///
			Position getProxyPort() const { return proxy_port_;}

			#ifdef BALL_QT_HAS_THREADS
			/// QWaitCondition to wake up threads, after Composites are unlocked
			QWaitCondition& getCompositesLockedWaitCondition() { return composites_locked_wait_condition_;}
			#endif

			/** Multithreaded code is used for serveral functions:
			    - Update of Representations
					- Simulations
					- Download PDB files
					To debug such code it is often usefull to to be able to run it in
					a singlethreaded mode. Every piece of multithreaded code should
					therefore call this method and decide if it should run without multiple threads.
					Furthermore most of the time, valid benchmark results can only be achived 
					with one single thread.
			*/
			bool useMultithreading()
				throw();

			/// See above
			void setMultithreading(bool state)
				throw() { multi_threading_mode_ = state;}

			//@}
			/**	@name	Debugging and Diagnostics
			*/
			//@{

			/** Internal state dump.
					Dump the current internal state of this mainControl to 
					the output ostream <b>s</b> with dumping depth <b>depth</b>.
					\param   s output stream where to output the internal state 
					\param   depth the dumping depth
			*/
			virtual void dump(std::ostream& s = std::cout, Size depth = 0) const
				throw();
					
			/** Open a file.
			 		To be derived from...
			*/
			virtual void openFile(const String& /*file*/) throw() {};

			///
			void saveBALLViewProjectFile(const String& filename);
			
			///
			void loadBALLViewProjectFile(const String& filename) throw();

			///
			void quickSave();

			///
			void quickLoad();

			///
			void processEvents(Size ms);

			//@}
			
			protected slots:

			/*_ This slot is called internally whenever the apply button
					of the Preferences dialog	is pressed.
					It calls among other things the method applyPreferences().
			*/
			virtual void applyPreferencesClicked_();

			//_ Called by timer to clear the text in the statusbar
			void clearStatusBarText_();

			// Connected to the delete entry
			virtual void deleteClicked();

			protected:

			virtual void initializePreferencesTab_()
				throw();

			//_  Called after receiving an SimulationThreadFinished event
			void stopedSimulation_();

			///
			void lockComposites_();

			/*_ Remove a composite.
					Every Representation, which was created for the Composite is deleted, by sending a 
					RepresentationMessage with type RepresentationMessage::REMOVE.\par
					Redraws representations of the parent of the Composite, if wished.
					\return bool <tt>true</tt> if the CompositeManager has the Composite
			*/
			bool remove_(Composite& composite, bool update_representations_of_parent = true, 
																				 bool to_delete = true)
				throw();

			/*_	Create a unique item ID for a menuentry by adding 1 to current_id_
			*/
			int getNextID_()
				throw();

			void selectRecursive_(Composite* composite)
				throw();

			void deselectRecursive_(Composite* composite)
				throw();

			/*_ Select the composite parents of the geometric objects.
					The GeometricObjectSelectionMessage is sent by the Scene.
			 */
			void selectComposites_(GeometricObjectSelectionMessage& message)
				throw();

			void reduceSelection_(Composite* const composite);

			//_ Called by constructors
			void setup_()
				throw();

			void complementSelectionHelper_(Composite& c);

			/** Show a busy cursor and a busy icon in the statusbar.
			*/
			void setBusyMode_(bool state);

			//_
			void setPreferencesEnabled_(bool state);

			void init_();

			//_
			FragmentDB fragment_db_;

			/*_ List with the selected composites
			*/
			HashSet<Composite*> 				selection_;

			/*_ List with the selected composites of the control.
					(Not the one with the checkboxes!)
			*/
			List<Composite*>						control_selection_;		

			/*_ Message label in the statusbar
					\see setStatusbarText
			*/
			QLabel* 										message_label_;

			PrimitiveManager 						primitive_manager_;
			CompositeManager 						composite_manager_;

			MainControlPreferences* 		main_control_preferences_;
			NetworkPreferences* 				network_preferences_;
			Preferences*								preferences_dialog_;
			int 			 									preferences_id_;
			int 			 									delete_id_;
			INIFile		 									preferences_file_;
			
			static int 									current_id_;
			bool 												composites_locked_;
			ModularWidget*							locking_widget_;
			bool 											  stop_simulation_;

			SimulationThread* 					simulation_thread_;

			bool 												multi_threading_mode_;

			/*_	A list containing all modular widgets.
					This list is modified by addModularWidget and
					removeModularWidget.
			*/
			List<ModularWidget*>				modular_widgets_;

			HashMap<Index, String>      menu_entries_hints_;

			QLabel*             simulation_icon_;
			static const char  *simulation_running_xpm_[];
			static const char  *simulation_stoped_xpm_[];

			String 							working_dir_;

			String 							logging_file_name_;
			bool 								logging_to_file_;
			File 								logging_file_;

			bool 								about_to_quit_;
			bool 								important_text_in_statusbar_;
			QTimer 							timer_;
			#ifdef 	BALL_QT_HAS_THREADS
			QMutex 							composites_locked_mutex_;
			QWaitCondition 			composites_locked_wait_condition_;
			#endif

			String 							proxy_;
			Position 						proxy_port_;

			Index 							stop_simulation_id_, complement_selection_id_, open_id_, save_project_id_;

			Index last_highlighted_menu_entry_;
};

#		ifndef BALL_NO_INLINE_FUNCTIONS
#			include <BALL/VIEW/KERNEL/mainControl.iC>
#		endif 
    
		}	// namespace VIEW
	} // namespace BALL

#endif // BALL_VIEW_KERNEL_MAINCONTROL_H
