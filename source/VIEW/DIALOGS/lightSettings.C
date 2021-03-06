// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// $Id: lightSettings.C,v 1.20.2.9 2005/12/10 13:54:40 amoll Exp $
//

#include <BALL/VIEW/DIALOGS/lightSettings.h>
#include <BALL/VIEW/WIDGETS/scene.h>

#include <qpushbutton.h>
#include <qlineedit.h> 
#include <qlistbox.h>
#include <qlabel.h>
#include <qbuttongroup.h>
#include <qradiobutton.h>
#include <qslider.h>

namespace BALL
{
	namespace VIEW
	{

LightSettings::LightSettings(QWidget* parent, const char* name, WFlags fl)
  : LightSettingsData(parent, name, fl),
		PreferencesEntry(),
		ignore_(false),
		current_light_(-1)
{
	if (parent == 0 || !RTTI::isKindOf<Scene>(*parent)) 
	{
		Log.error() << "LightSettings dialog must be created with a Scene as parent!" << std::endl;
		return;
	}

	stage_ = (dynamic_cast<Scene*>(parent))->getStage();
	if (stage_ == 0) 
	{
		Log.error() << "LightSettings dialog was created with a Scene as parent, which has no Stage!" << std::endl;
		return;
	}
	
	relative_to_camera->setChecked(true);
	defaultsPressed();
	setWidgetStackName("Lighting");
	registerWidgetForHelpSystem_(this, "scene.html#lightsources");
}


void LightSettings::updateFromStage()
	throw()
{
	lights_.clear();
	List<LightSource>::ConstIterator it = stage_->getLightSources().begin();
	for (; it != stage_->getLightSources().end(); it++)
	{
		lights_.push_back(*it);
	}

	update();
}


void LightSettings::update()
	throw()
{
	if (!lights_.size()) 
	{
		clearFields_();
		return;
	}

	if (lights_.size() != lights_list->count())
	{
		clearFields_();
	
		for (Position light_nr = 0; light_nr < lights_.size(); light_nr++)
		{
			lights_list->insertItem((String("Light ") + String(light_nr + 1)).c_str());
		}
	}

	if (getCurrentLightNumber_() == -1)
	{
		ignore_ = true;
		lights_list->setSelected(lights_.size() - 1, true);
		ignore_ = false;
		return;
	}
	
	getValues_();
	ignore_ = false;
}


void LightSettings::addLightPressed()
{
	if (lights_.size() >= 8) return;

	saveSettingsToLight_();

	LightSource light;
	light.setIntensity((float) 0.8);
	light.setType(LightSource::POSITIONAL);

	// create a kind of headlight
	// position light 20 space units behind camera position
	
	light.setRelativeToCamera(true);
	light.setPosition(Vector3(0, 4, -20));
	light.setDirection(Vector3(0, 0, 1));

	lights_.push_back(light);
	update();
}


void LightSettings::colorPressed()
{
	chooseColor(color_sample);
}


void LightSettings::defaultsPressed()
{
	lights_list->clear();
	lights_.clear();
	addLightPressed();
}


// store settings for last selected light
void LightSettings::saveSettingsToLight_()
	throw()
{
	if (current_light_ == -1 ||
			current_light_ >= (Index)lights_.size()) 
	{
		return;
	}

	LightSource& light = lights_[current_light_];
	light.setColor(color_sample->backgroundColor());

	try
	{
		Vector3 pos = getPosition_();
		Vector3 dir = getDirection_();
		bool relative = LightSettingsData::relative->isChecked();
		light.setRelativeToCamera(relative);

		// position and direction
		
		if (!relative)
		{
			Vector3 diff = dir - pos;
			dir = diff;
		}
		
		light.setPosition(pos);
		light.setDirection(dir);

		light.setIntensity((float)(intensity->value()) / 100.0);

		/////////////////////////////////////////////////////
		// type of light

		if 			(light_type->selected() == ambient)
		{
			light.setType(LightSource::AMBIENT);
		}
		else if (light_type->selected() == point)
		{
			light.setType(LightSource::POSITIONAL);
		}
		else
		{
			light.setType(LightSource::DIRECTIONAL);
		}
	}
	catch (Exception::GeneralException e)
	{
		Log.error() << "Invalid values in LightSettingsDialog" << std::endl;
		Log.error() << e;
	}
}
	

void LightSettings::lightSelected()
{
	if (!ignore_) saveSettingsToLight_();
	current_light_ = getCurrentLightNumber_();
	getValues_();
}


void LightSettings::removeLightPressed()
{
	Index current = getCurrentLightNumber_();
	if (current == -1) return;

	vector<LightSource>::iterator it = lights_.begin();
	for (Index i = 0; it != lights_.end() && i < current; it++)
	{
		i++;
	}
	lights_.erase(it);
	update();
}


void LightSettings::typeSelected()
{
	typeSelected_(light_type->selectedId());
}

void LightSettings::typeSelected_(Position type)
{
	bool is_ambient = (type == LightSource::AMBIENT);
	
	bool pos_enabled = type != LightSource::DIRECTIONAL && !is_ambient;

	position_x->setEnabled(pos_enabled);
	position_y->setEnabled(pos_enabled);
	position_z->setEnabled(pos_enabled);

	direction_x->setEnabled(!is_ambient);
	direction_y->setEnabled(!is_ambient);
	direction_z->setEnabled(!is_ambient);
	relative_to_camera->setEnabled(!is_ambient);
	not_relative->setEnabled(!is_ambient);
}

void LightSettings::getValues_()
	throw()
{
	Index current = getCurrentLightNumber_();
	if (current == -1) return;

	setControlsEnabled_(true);

	LightSource& light = lights_[current];

	color_sample->setBackgroundColor(light.getColor().getQColor());

	Vector3 pos = light.getPosition();
	Vector3 dir = light.getDirection();

	if (light.isRelativeToCamera())
	{
		relative->setChecked(true);
	}
	else
	{
		not_relative->setChecked(true);
		dir = pos + dir;
	}

	setPosition_(pos);
	setDirection_(dir);

	typeSelected_(light.getType());
	
	light_type->setButton(light.getType());
	intensity->setValue((Index)(light.getIntensity() * 100.0));
}


void LightSettings::clearFields_()
	throw()
{
	lights_list->clear();
	position_x->clear();
	position_y->clear();
	position_z->clear();
	direction_x->clear();
	direction_y->clear();
	direction_z->clear();
	setControlsEnabled_(false);
}


void LightSettings::setControlsEnabled_(bool state)
{
	remove_lights_button->setEnabled(state);
	position_x->setEnabled(state);
	position_y->setEnabled(state);
	position_z->setEnabled(state);
	direction_x->setEnabled(state);
	direction_y->setEnabled(state);
	direction_z->setEnabled(state);
	intensity->setEnabled(state);
	intensity->setValue(0);
	remove_lights_button->setEnabled(state);
	relative_to_camera->setEnabled(state);
	not_relative->setEnabled(state);
	color_button->setEnabled(state);
	light_type->setEnabled(state);
}

void LightSettings::apply()
	throw()
{
	saveSettingsToLight_();
	stage_->clearLightSources();
	for (Position p = 0; p < lights_.size(); p++)
	{
		stage_->addLightSource(lights_[p]);
	}
}


void LightSettings::intensityChanged()
{
	intensity_label->setText(String(intensity->value()).c_str());
}


void LightSettings::restoreDefaultValues(bool /*all*/)
	throw()
{
	defaultsPressed();
	color_sample->setBackgroundColor(white);
	lights_list->setCurrentItem(0);
}

void LightSettings::positionTypeChanged()
{
	if (getCurrentLightNumber_() == -1 || ignore_) 
	{
		return;
	}

	try
	{
		Vector3 pos = getPosition_();
		Vector3 dir = getDirection_();

		const Vector3& vp = stage_->getCamera().getViewPoint();

		if (relative->isChecked())
		{
			Vector3 diff = dir - pos;
			dir = stage_->calculateRelativeCoordinates(diff);

			pos -= vp;
			pos = stage_->calculateRelativeCoordinates(pos);
		}
		else
		{
			pos = stage_->calculateAbsoluteCoordinates(pos) + vp;
			dir = stage_->calculateAbsoluteCoordinates(dir) + pos;
		}

		setPosition_(pos);
		setDirection_(dir);
	}
	catch(...)
	{
		BALLVIEW_DEBUG
	}
}

void LightSettings::setPosition_(const Vector3& v)
{
	position_x->setText(createFloatString(v.x, 2).c_str());
	position_y->setText(createFloatString(v.y, 2).c_str());
	position_z->setText(createFloatString(v.z, 2).c_str());
}

void LightSettings::setDirection_(const Vector3& v)
{
	direction_x->setText(createFloatString(v.x, 2).c_str());
	direction_y->setText(createFloatString(v.y, 2).c_str());
	direction_z->setText(createFloatString(v.z, 2).c_str());
}

Vector3 LightSettings::getPosition_() 
	throw(Exception::InvalidFormat)
{
	return Vector3(String(position_x->text().ascii()).toFloat(),
				 			   String(position_y->text().ascii()).toFloat(),
			  				 String(position_z->text().ascii()).toFloat());
}

Vector3 LightSettings::getDirection_() 
	throw(Exception::InvalidFormat)
{
	return Vector3(String(direction_x->text().ascii()).toFloat(),
				 			   String(direction_y->text().ascii()).toFloat(),
			  				 String(direction_z->text().ascii()).toFloat());
}

Index LightSettings::getCurrentLightNumber_() const
{
	return lights_list->currentItem();
}

void LightSettings::restoreValues(bool)
{
	updateFromStage();
}
		
} } // NAMESPACE
