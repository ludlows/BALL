// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// $Id: displayProperties.sip,v 1.5.4.3 2005/12/02 14:01:22 amoll Exp $
//

class DisplayProperties
	: ModularWidget
{
%TypeHeaderCode
	#include <BALL/VIEW/DIALOGS/displayProperties.h>
	using namespace BALL;
	using namespace BALL::VIEW;
%End
  public:

// methods implemented by the BALL_EMBEDDED macro
static Size countInstances();
static DisplayProperties* getInstance(Position /index/);

//  DisplayProperties(QWidget *parent, const char*);
//	DisplayProperties(const DisplayProperties&)
~DisplayProperties() throw();
//  void onNotify(Message *message) throw();
//  void initializeWidget(MainControl&) throw();
//  void initializePreferencesTab(Preferences &preferences) throw();
//  void finalizePreferencesTab(Preferences &preferences) throw();
//  void applyPreferences() throw();
//  void checkMenu(MainControl&) throw();
  void createRepresentationMode();
  void modifyRepresentationMode(Representation*);
	Representation* createRepresentation(const PyCompositeList);
	bool getSettingsFromString(const String&) throw();
	void enableCreationForNewMolecules(bool) throw();
//  void show();
  void selectModel(int);
  void selectMode(int);
  void selectColoringMethod(int);
	void setSurfaceDrawingPrecision(float);
	void setDrawingPrecision(int);
	void setTransparency(int);
	void setCustomColor(const ColorRGBA&);
  void apply();
  void editColor();
  void editSelectionColor();
  void coloringOptionsPressed();
  void modelOptionsPressed();
  void precisionBoxChanged(int);
//  void transparencySliderChanged();
  void precisionSliderChanged();
	void coloringUpdatesChanged();
	void modelUpdatesChanged();
};
