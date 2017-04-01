// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// $Id: preferences.h,v 1.16.6.3 2005/09/29 14:01:31 amoll Exp $
//

#ifndef BALL_VIEW_DIALOGS_PREFERENCES_H
#define BALL_VIEW_DIALOGS_PREFERENCES_H

#include <qtabdialog.h>
#include <BALL/DATATYPE/string.h>
#include <BALL/VIEW/UIC/preferencesData.h>

#ifndef BALL_DATATYPE_HASHSET_H
# include <BALL/DATATYPE/hashSet.h>
#endif

#ifndef BALL_DATATYPE_HASHMAP_H
# include <BALL/DATATYPE/hashMap.h>
#endif

namespace BALL
{
	class INIFile;

	namespace VIEW
	{
		class PreferencesEntry;

		/** Preferences is a dialog for other preference pages to be
				inserted. An interface for adding and removing the pages is
				available. The Preferences dialog will be created by the MainControl. 
				All ModularWidget classes can add their own preferences pages.
				\ingroup ViewDialogs
		*/
		class BALL_VIEW_EXPORT Preferences 
			: public PreferencesData
		{
			Q_OBJECT
			
			public:

			/**	@name	Constructors
			*/	
			//@{

			/** Default Constructor.
					There are three buttons:
						-  apply - apply the preferences
						-  cancel - hide the preferences dialog
						-  defaults - reset the values of the currently selected page
					\param       parent the parent QWidget (See QTabDialog in the QT documentation)
					\param       name (See QDialog in the QT documentation)
			*/
			Preferences(QWidget *parent = NULL,	const char *name = NULL)
				throw();

			/// Copy constructur
			Preferences(const Preferences& preferences)
				throw();
			
			/// Destructor.
			virtual ~Preferences() throw();

			//@}
			/**	@name	Accessors
			*/
			//@{

			/** Check if pages are available.				
					\return bool <tt>true</tt> if pages are available
			*/
			bool hasPages()
				throw();

			/** Insert a new page.
					This method can be called inside ModularWidget::initializePreferencesTab
					\param  child a pointer to the new dialog
					\param  name the name of the new dialog
					\see    removePage
			*/
			void insertEntry(PreferencesEntry *child)
				throw();

			/** Remove a previously inserted page.
					This method can called inside ModularWidget::finalizePreferencesTab 
					\param  child a pointer to the dialog to be removed 
					\see    insertPage
			*/
			void removeEntry(PreferencesEntry *child)
				throw();

			/** Fetch the preferences (the position) from the INIFile <tt>inifile</tt>.
					This method will be called inside MainControl::show().
					\param  inifile the INIFile that contains the needed information
					\see    writePreferences
			*/
			void fetchPreferences(INIFile &inifile)
				throw();

			/** Write the preferences (the position) to the INIFile <tt> inifile</tt>.
					This method will be called in MainControl::aboutToExit.
					\param  inifile the INIFile to be written into
					\see    fetchPreferences
			*/
			void writePreferences(INIFile &inifile)
				throw();

			///
			const QWidget* currentEntry() const
				throw();

			///
			const QWidget* currentPage() const
				throw();

			public slots:

			//@}
			/** @name Public slots
			*/
			//@{

			/** Open the preferences dialog.
			*/
			void show();
			
			/// Show the given page
			void showEntry(QWidget* page);
						
			///
			void entrySelected(QListViewItem* item);

			///
			void setDefaultValues();

			///
			void cancelPreferences();

			/// Store the last applied values for undo
			void applyPreferences();

			///
			void showHelp();

			protected:

			void removeItem_(QListViewItem* item, bool update)
				throw();

			HashSet<PreferencesEntry*> entries_;
			HashMap<QListViewItem*, QWidget*> item_to_widget_;
			HashMap<QListViewItem*, PreferencesEntry*> item_to_entry_;
			HashMap<QWidget*, QListViewItem*> widget_to_item_;

			//@}
		};

  } // namespace VIEW
} // namespace BALL

#endif // BALL_VIEW_DIALOGS_PREFERENCES_H
