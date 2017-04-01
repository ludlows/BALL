// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// $Id: displayProperties.h,v 1.39.2.8 2005/12/02 14:01:25 amoll Exp $
//

#ifndef BALL_VIEW_DIALOGS_DISPLAYPROPERTIES_H
#define BALL_VIEW_DIALOGS_DISPLAYPROPERTIES_H

#ifndef BALL_VIEW_KERNEL_MODULARWIDGET_H
# include <BALL/VIEW/KERNEL/modularWidget.h>
#endif

#ifndef BALL_VIEW_DATATYPE_COLORRGBA_H
# include <BALL/VIEW/DATATYPE/colorRGBA.h>
#endif

#ifndef BALL_VIEW_KERNEL_PREFERENCESENTRY
# include <BALL/VIEW/KERNEL/preferencesEntry.h>
#endif

#include <BALL/VIEW/UIC/displayPropertiesData.h>

#include <qtimer.h>

namespace BALL
{
	class Composite;

	namespace VIEW
	{
		class Representation;
		
		class ColoringSettingsDialog;
		class ModelSettingsDialog;
		class Preferences;

		/**	Dialog for creating and changing representations for a selection of molecular objects.
				It can create a new Representation for a selection of Composite 's from the 
				MolecularControl. If a Representation is selected in the GeometricControl, it
				can be modified with this dialog.
				With the help of various combo boxes it is possible to customize the look of
				the graphical visualization (the model, the drawing precision, the drawing mode,
				the coloring method and the custom color).
				For information about the drawing precision see Representation. <br>
				If this dialog is used, it should be created with MainControl as parent. <br>
				If you want to add a new coloring method or type of Representation, have a look at
				VIEW/KERNEL/common.h.
				\ingroup ViewDialogs
		*/
		class BALL_VIEW_EXPORT DisplayProperties 
			: public DisplayPropertiesData,
				public ModularWidget,
				public PreferencesEntry
		{
			Q_OBJECT
				
			public:
			
			BALL_EMBEDDABLE(DisplayProperties,ModularWidget)

			/**	@name	Constructors and Destructors
			*/	
			//@{

			/** Default Constructor.
					Calls ModularWidget::registerWidget.
			*/
			DisplayProperties(QWidget *parent = NULL, const char* name = NULL)
				throw();

			/// Copy constructor just implemented for Python Interface, dont use it! 
			DisplayProperties(const DisplayProperties& dp)
			 throw();

			/** Destructor
			*/
			virtual ~DisplayProperties()
				throw();

			//@} 
			/**	@name	Accessors: inspectors and mutators 
			*/ 
			//@{

			/** Message handling method.
					Handles messages sent by other registered ConnectionObject objects.
					If a CompositeMessage with type NEW_MOLECULE is catched,
					the chosen graphical visualization
					will be applied to the Composite object and the follwing Message
					objects will be sent through the ConnectionObject tree:
						- CompositeMessage with type CENTER_CAMERA
						- RepresentationMessage with type NEW
					\par
					\param message the pointer to the message that should be processed
			*/
			virtual void onNotify(Message *message)
					throw();
			
			//@} 
			/**	ModularWidget methods 
			*/ 
			//@{
			
			/**	Initialize the popup menu <b>Display</b> with the entry
					<b>Display Properties</b>, which opens the dialog.
					This method is called automatically	immediately before the main application is started
					by MainControl::show()
					\param main_control the MainControl object to be initialized 
			*/
			virtual void initializeWidget(MainControl& main_control)
					throw();
				
			/// Insert the ModelSettingsDialog and the ColoringSettingsDialog into the Preferences
			virtual void initializePreferencesTab(Preferences &preferences)
				throw();

			/// Remove the ModelSettingsDialog and the ColoringSettingsDialog from the Preferences
			virtual void finalizePreferencesTab(Preferences &preferences)
				throw();

			///
			void applyPreferences()
				throw();

			/**	Menu checking method.
					This method is called by MainControl::checkMenus before a popup menu is shown.
					The menu entry <b>Display Properties</b> will be checked if this dialog is visible. 
					\param main_control the MainControl object whose menus should be checked
			*/
			virtual void checkMenu(MainControl& main_control)
					throw();

			/// Switch to the mode, that a new Representation will be created
			void createRepresentationMode();

			/// Switch to the mode, that an existing Representation will be modified
			void modifyRepresentationMode(Representation* rep);
	
			/// Settings from String
			bool getSettingsFromString(const String& data)
				throw();

			/// Set if Representations are automaticaly created for new Molecules
			void enableCreationForNewMolecules(bool state) 
				throw() { create_representations_for_new_molecules_ = state;}

			/// Get the Representation on which DisplayProperties is working on
			Representation* getRepresentation() 
				throw() { return rep_;}

			/* 	Create the new representation for the selection in the MolecularControl or for a given List of Composites.
					Called by onNotify() after receiving CompositeMessage::NEW_MOLECULE and by apply().
					To insert a new type of model, this is the only method in DisplayProperties you have to
					change (See also VIEW/KERNEL/common.h).
			*/
			virtual Representation* createRepresentation(const List<Composite*>& composites);
	
			public slots:
					
			//@} 
			/** @name Public slots 
			*/ 
			//@{

			/** Starts the dialog.
					Calls QDialog::raise().
			*/
			void show();

			
			/** Changes the model.
					This slot is connected to the model combo box and will be automatically
					called if the content of this combo box is changed.
					\param  index the position of the entry in the combobox
			*/
			void selectModel(int index);

			/** Changes the drawing mode.
					This slot is connected to the mode combo box and will be automatically
					called if the content of this combo box is changed.
					\param  index the position of the entry in the combobox
			*/
			void selectMode(int index);

			/** Changes the coloring method.
					This slot is connected to the coloring method combo box and will be automatically
					called if the content of this combo box is changed.
					\param  index the position of the entry in the combobox
			*/
			void selectColoringMethod(int index);

			///
			void setSurfaceDrawingPrecision(float value);
			
			///
			void setDrawingPrecision(int value);

			///
			void setTransparency(int value);

			///
			void setCustomColor(const ColorRGBA& color);
			
			/** Indicates the apply button was pressed.
					Applies the selected model with its selected properties to 
					the selected Composite objects or modifies an existing Representation.
					A SceneMessage will be sent to inform the Scene.
					A RepresentationMessage notifies the GeometricControl.
			*/
			virtual void apply();
			
			/** Opens the dialog for editing the custom color.
					Opens a QColorDialog from the QT-library.
			 */ 
			void editColor();

			/** Opens the color dialog for the color of selected items.
					\see BALL_SELECTED_COLOR
			*/
			void editSelectionColor();

			///
			void coloringOptionsPressed();

			///
			void modelOptionsPressed();

			///
			void precisionBoxChanged(int index);

			///
			void transparencySliderChanged();

			///
			void precisionSliderChanged();

			///
			void coloringUpdatesChanged();

			///
			void modelUpdatesChanged();

			//@}
			
			
			protected:
			
			//_ Set buttons and slider according to the values
			void checkDrawingPrecision_()
				throw();

			//_
			virtual void getAdvancedModelOptions_()
				throw();

			//_
			virtual void getAdvancedColoringOptions_()
				throw();

			//_
			virtual void applyModelSettings_(Representation& rep);
			
			//_
			virtual void applyColoringSettings_(Representation& rep);

			//_
			bool isNotBusy_();

			protected slots:
			//_
			void checkMenu_();

			protected:
			
			// --------------------------------------------------------------------------------
			// attributs
			// --------------------------------------------------------------------------------
			ModelSettingsDialog* 			model_settings_;
			ColoringSettingsDialog* 	coloring_settings_;
			Preferences* 							preferences_;
			
			// the menu entry id of the dialog
			int 						id_;
			
			// used by GeometricControl to modify an existing representation
			Representation* rep_;

			ColorRGBA 			custom_color_;
			bool 						advanced_options_modified_;
			bool 						create_representations_for_new_molecules_;
			bool 						changed_selection_color_;
			QTimer 					timer_;
		};

} } // namespaces

#endif // BALL_VIEW_DIALOGS_DISPLAYPROPERTIES_H
