// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// $Id: modularWidget.h,v 1.19.6.8 2005/09/29 14:01:33 amoll Exp $
//

#ifndef BALL_VIEW_WIDGETS_MODULARWIDGET_H
#define BALL_VIEW_WIDGETS_MODULARWIDGET_H

#ifndef BALL_CONCEPT_EMBEDDABLE_H
#	include <BALL/CONCEPT/embeddable.h>
#endif

#ifndef BALL_VIEW_KERNEL_CONNECTIONOBJECT_H
# include <BALL/VIEW/KERNEL/connectionObject.h>
#endif

class QObject;
class QMenuBar;
class QWidget;

namespace BALL
{
	class INIFile;
	class FragmentDB;

	namespace VIEW
	{
		class Preferences;
		class MainControl;

		/**	ModularWidget is a base class.
				BALLVIEW provides a simple mechanism for Rapid Application Prototyping based 
				on <b>modular widgets</b>.	Using the modular widgets and the MainControl class
				(or more precisely, classes derived from MainControl) applications can be developed
				in a plug&play style.	The most important parts of a visualization application are	
				implemented as classes derived from ModularWidget.
				Inserting an instance of one of these widgets automatically creates the required
				menus in the menubar of the main window, sets up the required connections,
				and registers the class instance for use from the embedded Python interpreter 
				(if desired).
				<br>
				When implementing classes derived from ModularWidget and you want to give access 
				to the class from Python, please remember	to include the BALL_EMBEDDABLE 
				macro in the public section	of your class declaration. Also make sure that the QT 
				class you	derive from (e.g. QWidget) is the <b>first</b> base class and 
				ModularWidget	second. 
				<br>
				<b>Remember:</b> A ModularWidget is not notified by the Messages it sends itself!
				\see MainControl
				\see Embeddable
				\see PyInterpreter
				\see PyWidget
			\ingroup ViewKernelConnectivity
		*/
		class BALL_VIEW_EXPORT ModularWidget
			: public Embeddable,	
				public ConnectionObject
		{
		  public:
			
			BALL_EMBEDDABLE(ModularWidget,Embeddable)
			BALL_CREATE(ModularWidget)
			
			/**	@name	Constructors
			*/	
			//@{

			/** Default Constructor.
					Set the name of this ModularWidget to <b>name</b>.
					\param      name the name of this modularWidget
			*/
			ModularWidget(const char* name = "<ModularWidget>")
				throw();
				
			/**	Copy constructor.
			*/
			ModularWidget(const ModularWidget& widget)
				throw();
				
			//@}
			/** @name Destructors 
			*/
			//@{

			/** Destructor.
			 		If the ModularWidget was registered, it will call MainControl::removeModularWidget.
			*/
			virtual ~ModularWidget()
				throw();

			///
			virtual void destroy()
				throw() { clear(); };

			/** Explicit default initialization.
			 		Currently does nothing.
			*/
			virtual void clear()
				throw();

			//@}	
			/**	@name	Accessors: inspectors and mutators 
			*/
			//@{
			
			/**	Register the widget <b>mwidget</b> to the MainControl. ModularWidget
					objects must always be created with MainControl as parent and must have
					this method in their constructors. This method connects them to the
					MainControl object.
					\param  mwidget the ModularWidget to be registered to the MainControl
			*/
			static void registerWidget(ModularWidget* mwidget)
				throw(Exception::NullPointer);
				
			/**	Initialize the widget.
					This method is called automatically	immediately before the main application 
					is started. It should add the widget's menu entries and connections (if required).
					This method will be called by MainControl::show.
					\param main_control the MainControl object to be initialized with this ModularWidget
					\see   finalizeWidget()
			*/
			virtual void initializeWidget(MainControl& main_control);
			
			/**	Remove the widget custom items, e.g all menu entries.
					This method should reverse all actions performed in initializeWidget
					(remove menu entries and connections of this ModularWidget).
					Call this method also in derived classes finalizeWidget to remove the menu entries.
					This method will be called by MainControl::aboutToExit().
					\param main_control the MainControl object to be finalized with this ModularWidget
					\see   initializeWidget
			*/
			virtual void finalizeWidget(MainControl& main_control);
			
			/**	Menu checking method.
					This method is called MainControl::checkMenus before a popup menu is shown.
					It should be used to update the state of menu entries (e.g. disable or enable entries).
					\param main_control the MainControl object whose menus should be checked
			*/
			virtual void checkMenu(MainControl& main_control)
				throw();
			
			/** Initialize a preferences tab for the widget (if needed).
					This method can be used to create preferences widgets that can be inserted
					into the Preferences dialog with the method insertTab.
					This method is called automatically by MainControl::show at the start of the application.
					\param preferences the Preferences dialog of the MainControl
			*/
			virtual void initializePreferencesTab(Preferences& preferences)
				throw();
			
			/**	Remove the preferences tab.
					This method can remove a preferences widget (if created in 
					initializePreferencesTab)	from the Preferences dialog of the MainControl.
					This method is called automatically by MainControl::aboutToExit() at the end of the application.
					\param  preferences the Preferences dialog of the MainControl
			*/
			virtual void finalizePreferencesTab(Preferences& preferences)
				throw();
			
			/** Apply the preferences of the specific tab.
					In this method the widget can extract any changed values from
					its preferences tab (if required).
					This method is called automatically by the applyPreferencesTab from the
					MainControl object if the apply button in the Preferences dialog
					is pressed.
					\param  preferences the Preferences dialog of the MainControl
					\see    initializePreferencesTab
					\see    finalizePreferencesTab
					\see    applyPreferencesTab
			*/
			virtual void applyPreferences()
				throw() {};

			/** Fetch the widgets preferences from the INIFile.
					This method is called automatically by MainControl::show() at the start of the application.
					\param  inifile the INIFile that contains the needed values
					\see    writePreferences
			*/
			virtual void fetchPreferences(INIFile& inifile)
				throw();
			
			/** Writes the widgets preferences to the INIFile.
					This method is called automatically by MainControl::aboutToExit at the end of the application.
					\param  inifile the INIFile to contain the values
					\see    fetchPreferences
			*/
			virtual void writePreferences(INIFile& inifile)
				throw();

			/** Set the text of the statusbar of the main application.
			 		<b>Note:</b> The ModularWidget must be registered to a MainControl.
					Implemented for convenience.
			 */
			virtual void setStatusbarText(String text, bool important = false)
				throw();

			/** Return the MainControl of this ModularWidget
					Implemented for convenience.
			*/
			MainControl* getMainControl() const
				throw();

			///	Implemented for convenience.
			String getWorkingDir()
				throw();

			///	Implemented for convenience.
			void setWorkingDir(const String& dir)
				throw();

			/** Return the FragmentDB.
			 		<b>Note:</b> The ModularWidget must be registered to a MainControl.
					Implemented for convenience.
			*/
			FragmentDB& getFragmentDB() const
				throw();

			/** Try to get an exclusive lock on the Composites, so that they can not be altered by
			 		any other ModularWidget.
			*/
			bool lockComposites()
				throw();

			/// Unlock the Composites.
			bool unlockComposites()
				throw();

			/// Wrapper for MainControl::menuBar()
			QMenuBar* menuBar() 
				throw();

			Index insertMenuEntry (Index parent_id, const String& name, const QObject* receiver = 0, 
													 const char* slot = 0, Index accel = 0, Index pos = -1)
				throw();

			///
			void setMenuHint(const String& hint);

			///
			void setMenuHelp(const String& url);

			///
			virtual void registerWidgetForHelpSystem(const QWidget* widget, const String& url);

			///
			virtual void registerMenuEntryForHelpSystem(Index entry, const String& docu_entry);

			//@}
			/**	@name	Debugging and Diagnostics
			*/
			//@{

			/** Internal state dump.
					Dump the current internal state of this mainControl to 
					the output ostream <b>s</b> with dumping depth <b>depth</b>.
					Calls ConnectionObject::dump.
					\param   s output stream where to output the internal state 
					\param   depth the dumping depth
			*/
			virtual void dump(std::ostream& s = std::cout, Size depth = 0) const
				throw();

			//@}

			void setWorkingDirFromFilename_(String filename)
				throw();

			void removeMenuEntries();

			virtual void showHelp(const String& url);

			protected:

			//_ id in the menubar entry "WINDOWS" for every widget
			Index window_menu_entry_id_;

			//_ should there be an entry to switch the window on and off?
			bool show_window_enty_;

			//_ should the widget be visible, if no config file entry exists?
			bool default_visible_;

			vector<std::pair<Index, Index> > menu_ids_;

			Index last_parent_id_, last_id_;
		}; 
  
	} // namespace VIEW
} // namespace BALL

#endif // BALL_VIEW_WIDGETS_MODULARWIDGET_H
