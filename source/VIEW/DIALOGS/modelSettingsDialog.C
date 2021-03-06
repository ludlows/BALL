// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// $Id: modelSettingsDialog.C,v 1.34.4.5 2005/12/10 13:54:40 amoll Exp $
//

#include <BALL/VIEW/DIALOGS/modelSettingsDialog.h>


#include <BALL/VIEW/MODELS/backboneModel.h>
#include <BALL/VIEW/MODELS/cartoonModel.h>
#include <BALL/VIEW/MODELS/ballAndStickModel.h>
#include <BALL/VIEW/MODELS/lineModel.h>
#include <BALL/VIEW/MODELS/surfaceModel.h>
#include <BALL/VIEW/MODELS/vanDerWaalsModel.h>
#include <BALL/VIEW/MODELS/HBondModel.h>
#include <BALL/VIEW/MODELS/forceModel.h>
#include <BALL/VIEW/MODELS/standardColorProcessor.h>

#include <BALL/DATATYPE/string.h>
#include <BALL/FORMAT/INIFile.h>

#include <qslider.h>
#include <qlabel.h>
#include <qwidgetstack.h>
#include <qradiobutton.h>
#include <qbuttongroup.h>

namespace BALL
{
	namespace VIEW
	{

		ModelSettingsDialog::ModelSettingsDialog( QWidget* parent,  const char* name, WFlags fl )
			: ModelSettingsDialogData( parent, name, fl ),
				PreferencesEntry()
		{
			setINIFileSectionName("MODEL_OPTIONS");

			setDefaultValues_();

			registerObject_(stick_radius_slider);

			registerObject_(ball_stick_cylinder_radius_slider);
			registerObject_(ball_stick_sphere_radius_slider);
			registerObject_(ball_stick_dashed_bonds);
			
			registerObject_(vdw_radius_factor_slider);
			
			registerObject_(surface_probe_radius_slider);
			
			registerObject_(tube_radius_slider);
			
			registerObject_(cartoon_tube_radius_slider);
			registerObject_(cartoon_helix_radius_slider);
			registerObject_(strand_arrow_width_slider);
			registerObject_(strand_height_slider);
			registerObject_(strand_width_slider);
			registerObject_(cartoon_dna_helix_radius_slider);
			registerObject_(cartoon_dna_ladder_radius_slider);
			registerObject_(cartoon_dna_base_radius_slider);
			registerObject_(dna_cartoon_model_type);
			registerObject_(ribbons_enabled);
			registerObject_(two_colored_ribbons);
			
			registerObject_(force_scaling_slider);
			registerObject_(force_max_length_slider);

			registerObject_(hbonds_radius_slider);

			setWidgetStackName("Models");
			setWidgetStack(widget_stack);
		}

		void ModelSettingsDialog::setDefaultValues_()
		{
			AddBallAndStickModel dummy1;
			dummy1.enableStickModel();
			getSettings(dummy1);
			
			AddBallAndStickModel dummy2;
			getSettings(dummy2);
			
			AddVanDerWaalsModel dummy3;
			getSettings(dummy3);
			
			AddSurfaceModel dummy4;
			getSettings(dummy4);

			AddBackboneModel dummy5;
			getSettings(dummy5);

			AddCartoonModel dummy6;
			getSettings(dummy6);
		
			HBondModelProcessor dummy7;
			getSettings(dummy7);
			
			ForceModel dummy8;
			getSettings(dummy8);
		}

		float ModelSettingsDialog::getFloatValue_(const QSlider* const& slider) const
			throw()
		{
			return (float(slider->value())) / 10.0;
		}

		void ModelSettingsDialog::setValue_(QSlider* le, float value)
			throw()
		{
			le->setValue((Size)(value * 10.0));
		}

		void ModelSettingsDialog::setLabelText_(QLabel* label, const QSlider* const from)
			throw()
		{
			String s((float)from->value() / 10.0);
			s.trimRight("0");
			if (s.hasSuffix(".")) s+= "0";
			label->setText(s.c_str());
		}

		void ModelSettingsDialog::applySettingsTo(ModelProcessor& mp) const
			throw()
		{
			if (RTTI::isKindOf<AddLineModel>(mp))
			{
				return;
			}
					
			if (RTTI::isKindOf<AddBallAndStickModel>(mp))
			{
				AddBallAndStickModel& bsm = *((AddBallAndStickModel*)&mp);
				if (bsm.isStickModel())
				{
					bsm.setStickRadius(getStickStickRadius());
				}
				else
				{
					bsm.setStickRadius(getBallAndStickStickRadius());
					bsm.setBallRadius(getBallRadius());
				}
				bsm.enableDashedBonds(bsm.isBallAndStickModel() && ballAndStickDashedBondsEnabled());
				return;
			}
					
			if (RTTI::isKindOf<AddSurfaceModel>(mp))
			{
				((AddSurfaceModel*)&mp)->setProbeRadius(getSurfaceProbeRadius());
				return;
			}
					
			if (RTTI::isKindOf<AddVanDerWaalsModel>(mp))
			{
				((AddVanDerWaalsModel*) &mp)->setVDWRadiusFactor(getVDWRadiusFactor());
				return;
			}
					
			if (RTTI::isKindOf<AddCartoonModel>(mp))
			{
				AddCartoonModel& cm = *dynamic_cast<AddCartoonModel*>(&mp);
				cm.setTubeRadius(getCartoonTubeRadius());
				cm.setHelixRadius(getCartoonHelixRadius());
				cm.setArrowWidth(getCartoonArrowWidth());
				cm.setStrandHeight(getCartoonStrandHeight());
				cm.setStrandWidth(getCartoonStrandWidth());
				cm.setDrawDNAAsLadderModel(cartoon_dna_ladder->isChecked());
				cm.setDNALadderRadius(getDNALadderRadius());
				cm.setDNABaseRadius(getDNABaseRadius());
				cm.setDNAHelixRadius(getDNAHelixRadius());
				cm.enableRibbons(ribbons_enabled->isChecked());
				cm.enableTwoColors(two_colored_ribbons->isChecked());
				return;
			}

			// backbone model after cartoon model !!!
			if (RTTI::isKindOf<AddBackboneModel>(mp))
			{
				((AddBackboneModel*) &mp)->setTubeRadius(getTubeRadius());
				return;
			}
					
			if (RTTI::isKindOf<HBondModelProcessor>(mp))
			{
				((HBondModelProcessor*) &mp)->setRadius(getHBondsRadius());
				return;
			}
					
			if (RTTI::isKindOf<ForceModel>(mp))
			{
				((ForceModel*) &mp)->setMaxLength((float)(getForceMaxLength()) / 10.0);
				((ForceModel*) &mp)->setScaling((float)(getForceScaling()) / 10.0);
				return;
			}
		}


		ModelProcessor* ModelSettingsDialog::createModelProcessor(ModelType type) const
			throw(Exception::InvalidOption)
		{
			ModelProcessor* model_processor = 0;

			switch (type)
			{
				case MODEL_LINES:
					model_processor = new AddLineModel;
					break;
					
				case MODEL_STICK:
					model_processor = new AddBallAndStickModel;
					((AddBallAndStickModel*)model_processor)->enableStickModel();
					break;
					
				case MODEL_BALL_AND_STICK:
					model_processor = new AddBallAndStickModel;
					((AddBallAndStickModel*)model_processor)->enableBallAndStickModel();
					break;
					
				case MODEL_SE_SURFACE:
					model_processor = new AddSurfaceModel;
					((AddSurfaceModel*)model_processor)->setType(SurfaceProcessor::SOLVENT_EXCLUDED_SURFACE);	
					break;
					
				case MODEL_SA_SURFACE:
					model_processor = new AddSurfaceModel;
					((AddSurfaceModel*)model_processor)->setType(SurfaceProcessor::SOLVENT_ACCESSIBLE_SURFACE);	
					break;
					
				case MODEL_VDW:
					model_processor = new AddVanDerWaalsModel;
					((AddVanDerWaalsModel*) model_processor)->setVDWRadiusFactor(getVDWRadiusFactor());
					break;

				case MODEL_BACKBONE:
					model_processor = new AddBackboneModel;
					((AddBackboneModel*) model_processor)->setTubeRadius(getTubeRadius());
					break;

				case MODEL_CARTOON:
					model_processor = new AddCartoonModel;
					break;
					
				case MODEL_HBONDS:
					model_processor = new HBondModelProcessor;
					break;

				case MODEL_FORCES:
					model_processor = new ForceModel;
					break;
					
				default:
					throw(Exception::InvalidOption(__FILE__, __LINE__, type));
			}

			applySettingsTo(*model_processor);

			return model_processor;
		}


		void ModelSettingsDialog::getSettings(const ModelProcessor& mp)
			throw()
		{
			if (RTTI::isKindOf<AddBallAndStickModel>(mp))
			{
				if (((AddBallAndStickModel*) &mp)->isStickModel())
				{
					setStickStickRadius(((AddBallAndStickModel*)&mp)->getStickRadius());
				}
				else
				{
					setBallAndStickStickRadius(((AddBallAndStickModel*)&mp)->getStickRadius());
					setBallRadius(((AddBallAndStickModel*)&mp)->getBallRadius());
				}
				return;
			}
					
			if (RTTI::isKindOf<AddSurfaceModel>(mp))
			{
				setSurfaceProbeRadius(((AddSurfaceModel*)&mp)->getProbeRadius());
				return;
			}
					
			if (RTTI::isKindOf<AddVanDerWaalsModel>(mp))
			{
				setVDWRadiusFactor(((AddVanDerWaalsModel*) &mp)->getVDWRadiusFactor());
				return;
			}
					
			if (RTTI::isKindOf<AddCartoonModel>(mp))
			{
				AddCartoonModel& cm = *(AddCartoonModel*)(&mp);
				setCartoonTubeRadius(cm.getTubeRadius());
				setCartoonHelixRadius(cm.getHelixRadius());
				setCartoonArrowWidth(cm.getArrowWidth());
				setCartoonStrandHeight(cm.getStrandHeight());
				setCartoonStrandWidth(cm.getStrandWidth());
				setCartoonDNAHelixRadius(cm.getDNAHelixRadius());
				setCartoonDNALadderRadius(cm.getDNALadderRadius());
				setCartoonDNABaseRadius(cm.getDNABaseRadius());
				cartoon_dna_ladder->setChecked(cm.drawDNAAsLadderModel());
				cartoon_dna_wac->setChecked(!cm.drawDNAAsLadderModel());
				two_colored_ribbons->setChecked(cm.twoColorsEnabled());
				ribbons_enabled->setChecked(cm.ribbonsEnabled());
				return;
			}

			// after derived class
			if (RTTI::isKindOf<AddBackboneModel>(mp))
			{
				setTubeRadius(((AddBackboneModel*) &mp)->getTubeRadius());
				return;
			}
					
			if (RTTI::isKindOf<HBondModelProcessor>(mp))
			{
				setHBondRadius(((HBondModelProcessor*) &mp)->getRadius());
				return;
			}
					
			if (RTTI::isKindOf<ForceModel>(mp))
			{
				setForceScaling(((ForceModel*) &mp)->getScaling());
				setForceMaxLenght(((ForceModel*) &mp)->getMaxLength());
				return;
			}
		}

		QWidget* ModelSettingsDialog::getEntryFor(ModelType type)
			throw()
		{
			switch (type)
			{
				case MODEL_LINES:
					return 0;
					
				case MODEL_STICK:
					return widget_stack->widget(0);

				case MODEL_BALL_AND_STICK:
					return widget_stack->widget(1);
					
				case MODEL_SE_SURFACE:
				case MODEL_SA_SURFACE:
					return widget_stack->widget(3);
					
				case MODEL_VDW:
					return widget_stack->widget(2);

				case MODEL_BACKBONE:
					return widget_stack->widget(4);

				case MODEL_CARTOON:
					return widget_stack->widget(5);
					
				case MODEL_HBONDS:
					return widget_stack->widget(6);

				case MODEL_FORCES:
					return widget_stack->widget(7);

				default:
					break;
			}

			return 0;
		}


	} // namespace VIEW
} // namespace BALL
