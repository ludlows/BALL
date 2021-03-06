// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// $Id: coloringSettingsDialog.C,v 1.37.2.5 2005/12/10 13:54:39 amoll Exp $
//

#include <BALL/VIEW/DIALOGS/coloringSettingsDialog.h>
#include <BALL/VIEW/MODELS/standardColorProcessor.h>
#include <BALL/KERNEL/PTE.h>

#include <qcolordialog.h>
#include <qslider.h>
#include <qlabel.h>
#include <qwidgetstack.h>
#include <qcheckbox.h>

namespace BALL
{
	namespace VIEW
	{

		QColorTableItem::QColorTableItem(QTable *t, EditType et, const ColorRGBA& color)
			: QTableItem(t,et,""),
				color_rgba_(color)
		{ 
		}

		void QColorTableItem::paint(QPainter *p, const QColorGroup &cg, const QRect &cr, bool selected)
		{
			QColorGroup g( cg );
			g.setColor( QColorGroup::Base, QColor (color_rgba_.getQColor()));
			QTableItem::paint( p, g, cr, selected );
		}

		// ==============================================================================================
		QColorTable::QColorTable(QWidget* parent, const char* name) 
			throw()
			: QTable(parent, name),
				setting_content_(false)
		{
			setNumCols(2);
			horizontalHeader()->setLabel(1, "Color");
			setGeometry(5,5, 400, 388);
			setColumnWidth(1, 230);
			setSelectionMode(NoSelection);
		}

		void QColorTable::setNamesTitle(const String& name)
			throw()
		{
			horizontalHeader()->setLabel(0, name.c_str());
		}

		String QColorTable::getNamesTitle() const
			throw()
		{
			return horizontalHeader()->label(0).ascii();
		}
			
		void QColorTable::setContent(const vector<String>& names, const vector<ColorRGBA>& colors)
			throw()
		{
			setting_content_ = true;
			colors_ = colors;
			names_ = names;

			setNumRows(colors_.size());
			for (Position p = 0; p < names_.size(); p++)
			{
				QColorTableItem* c2 = new QColorTableItem(this, QTableItem::WhenCurrent, colors_[p]);
				setText(p, 0, names_[p].c_str());
				setItem(p, 1, c2 );
			}
			setting_content_ = false;
		}

		void QColorTable::setColors(const vector<ColorRGBA>& colors)
			throw()
		{
			setting_content_ = true;
			colors_ = colors;

			setNumRows(colors_.size());
			for (Position p = 0; p < names_.size(); p++)
			{
				QColorTableItem* c2 = new QColorTableItem(this, QTableItem::WhenCurrent, colors_[p]);
				setItem(p, 1, c2 );
			}
			setting_content_ = false;
		}

		QWidget* QColorTable::beginEdit(int row, int col, bool)
		{
			if (col == 0 || setting_content_) return 0;
			ColorRGBA old_rgba(((QColorTableItem*)item(row,col))->getColor());
			QColor qcolor = QColorDialog::getColor(old_rgba.getQColor());
			if (!qcolor.isValid()) return 0;

			ColorRGBA new_color(qcolor);
			((QColorTableItem*)item(row,col))->setColor(new_color);
			updateCell(row, col);
			colors_[row] = new_color;
			return 0;
		}
		
		bool QColorTable::getValue(String& value) const
		{
			value.clear();

			vector<ColorRGBA> colors = getColors();

			for (Position p = 0; p < colors.size(); p ++)
			{
				value += colors[p];
				value += ";";
			}

			return true;
		}

		bool QColorTable::setValue(const String& value)
		{
			vector<String> fields;
			vector<ColorRGBA> colors;

			Size nr = value.split(fields, ";");

			ColorRGBA color;

			for (Position p = 0; p < nr; p ++)
			{
				color = fields[p];
				colors.push_back(color);
			}
			
			if (colors.size() != getNames().size())
			{
				BALLVIEW_DEBUG;
				return false;
			}

			setColors(colors);

			return true;
		}

					
		// =========================================================================================
		ColoringSettingsDialog::ColoringSettingsDialog( QWidget* parent,  const char* name, WFlags fl )
			: ColoringSettingsDialogData(parent, name, fl),
				PreferencesEntry()
		{
			setINIFileSectionName("COLORING_OPTIONS");

			element_table_ = new QColorTable(widget_stack->widget(0), "Elements");
			residue_table_ = new QColorTable(widget_stack->widget(2), "Residues");
			chain_table_   = new QColorTable(widget_stack->widget(10), "Chains");
			molecule_table_= new QColorTable(widget_stack->widget(11), "Molecules");
			setDefaultValues_();

			registerObject_(first_residue_label);
			registerObject_(middle_residue_label);
			registerObject_(last_residue_label);
			
			registerObject_(negative_charge_label);
			registerObject_(neutral_charge_label);
			registerObject_(positive_charge_label);

			registerObject_(null_distance_label);
			registerObject_(max_distance_label);
			registerObject_(max_distance_slider);
			registerObject_(distance_show_selected);

			registerObject_(minimum_tf_label);
			registerObject_(maximum_tf_label);
			registerObject_(unassigned_tf_label);
			registerObject_(max_tf_slider);

			registerObject_(minimum_o_label);
			registerObject_(unassigned_o_label);
			registerObject_(maximum_o_label);

			registerObject_(helix_color_label);
			registerObject_(coil_color_label);
			registerObject_(strand_color_label);
			registerObject_(turn_color_label);

			registerObject_(force_max_color_label);
			registerObject_(force_min_color_label);
			registerObject_(force_max_value_slider);
			registerObject_(force_min_value_slider);

			registerObject_(acidic_color_label);
			registerObject_(aromatic_color_label);
			registerObject_(basic_color_label);
			registerObject_(hydrophobic_color_label);
			registerObject_(other_color_label);
			registerObject_(polar_color_label);

			registerObject_(element_table_);
			registerObject_(chain_table_);
			registerObject_(residue_table_);
			registerObject_(molecule_table_);

			setWidgetStackName("Model Colors");
			setWidgetStack(widget_stack);
		}

		void ColoringSettingsDialog::setDefaultValues_()
		{
			vector<String> 		names;
			vector<ColorRGBA> colors;

			// =============================================================
			// setting element colors
			{
				// create a dummy processor to get the default values
				ElementColorProcessor elp;
				const HashMap<Position, ColorRGBA>& color_hash_map = elp.getColorMap();
				HashMap<Position, ColorRGBA>::ConstIterator it = color_hash_map.begin();

				for(; it != color_hash_map.end(); it++)
				{
					if (it->first == 0) continue;
					names.push_back(PTE[it->first].getSymbol());
					colors.push_back(it->second);
				}

				names.push_back(PTE[0].getSymbol());
				colors.push_back(color_hash_map[0]);
				element_table_->setNamesTitle("Element");
				element_table_->setContent(names, colors);
				names.clear();
				colors.clear();
			}

			// =============================================================
			// setting residue name colors
			// create a dummy processor to get the default values
			{
				ResidueNameColorProcessor rcp;
				const StringHashMap<ColorRGBA>& color_map = rcp.getColorMap();
				StringHashMap<ColorRGBA>::ConstIterator it2 = color_map.begin();

				for(; it2 != color_map.end(); it2++)
				{
					names.push_back(it2->first);
					colors.push_back(it2->second);
				}

				residue_table_->setNamesTitle("Residue");
				residue_table_->setContent(names, colors);
			}
				
			// =============================================================
			getSettings(ResidueNumberColorProcessor());
			getSettings(AtomChargeColorProcessor());
			getSettings(AtomDistanceColorProcessor());
			getSettings(TemperatureFactorColorProcessor());
			getSettings(OccupancyColorProcessor());
			getSettings(SecondaryStructureColorProcessor());
			getSettings(ForceColorProcessor());
			getSettings(ResidueTypeColorProcessor());
			getSettings(ChainColorProcessor());
			getSettings(MoleculeColorProcessor());
		}

		vector<ColorRGBA> ColoringSettingsDialog::getColors(ColoringMethod method) const
			throw()
		{
			vector<ColorRGBA> colors;
			QColorTable* table = 0;
			switch (method)
			{
				case COLORING_ELEMENT:
				{
					table = element_table_; 
					if (table->numRows() > 0)
					{
						colors.push_back(((QColorTableItem*)table->item(table->numRows() - 1, 1))->getColor());
					}
					break;
				}
				case COLORING_RESIDUE_INDEX: 	table = residue_table_; break;
				case COLORING_CHAIN: 					table = chain_table_  ; break;
				case COLORING_MOLECULE: 			table = molecule_table_; break;
				default: return colors;
			}
					
			for (Position p = 0; p < (Position)table->numRows(); p++)
			{
				colors.push_back(((QColorTableItem*)table->item(p, 1))->getColor());
			}

			return colors;
		}

		void ColoringSettingsDialog::applySettingsTo(ColorProcessor& cp) const
			throw()
		{
			if (RTTI::isKindOf<CustomColorProcessor>(cp)) return;

			if (RTTI::isKindOf<ElementColorProcessor>(cp))
			{
				vector<ColorRGBA> colors = getColors(COLORING_ELEMENT);
				for (Position p = 0; p < colors.size(); p++)
				{
					(*(ElementColorProcessor*)&cp).getColorMap()[p] = colors[p];
				}
				return;
			}
			
			if (RTTI::isKindOf<ResidueNameColorProcessor>(cp))
			{
				for (Position p = 0; p < (Position)residue_table_->numRows(); p++)
				{
					(*(ResidueNameColorProcessor*)&cp).getColorMap()[residue_table_->item(p,0)->text().ascii()] = 
						((QColorTableItem*)residue_table_->item(p,1))->getColor();
				}
				return;
			}

			if (RTTI::isKindOf<ResidueNumberColorProcessor>(cp))
			{
				ResidueNumberColorProcessor& dp = (*(ResidueNumberColorProcessor*)&cp);
				dp.setFirstColor(getLabelColor_(first_residue_label));
				dp.setMiddleColor(getLabelColor_(middle_residue_label));
				dp.setLastColor(getLabelColor_(last_residue_label));
				return;
			}

			if (RTTI::isKindOf<AtomChargeColorProcessor>(cp))
			{
				AtomChargeColorProcessor& dp = (*(AtomChargeColorProcessor*)&cp);
				dp.getColors()[0] = (getLabelColor_(negative_charge_label));
				dp.getColors()[1] = (getLabelColor_(neutral_charge_label));
				dp.getColors()[2] = (getLabelColor_(positive_charge_label));
				return;
			}

			if (RTTI::isKindOf<AtomDistanceColorProcessor>(cp))
			{
				AtomDistanceColorProcessor& dp = (*(AtomDistanceColorProcessor*)&cp);
				dp.setNullDistanceColor(getLabelColor_(null_distance_label));
				dp.setMaxDistanceColor(getLabelColor_(max_distance_label));
				dp.setDistance(((float)max_distance_slider->value()) / 10.0);
 				dp.setShowSelected(distance_show_selected->isChecked());
				return;
			}

			if (RTTI::isKindOf<OccupancyColorProcessor>(cp))
			{
				OccupancyColorProcessor& dp = (*(OccupancyColorProcessor*)&cp);
				dp.getColors()[0] = (getLabelColor_(minimum_o_label));
				dp.getColors()[1] = (getLabelColor_(maximum_o_label));
				return;
			}

			if (RTTI::isKindOf<SecondaryStructureColorProcessor>(cp))
			{
				SecondaryStructureColorProcessor& dp = (*(SecondaryStructureColorProcessor*)&cp);

				dp.setHelixColor(getLabelColor_(helix_color_label));
				dp.setCoilColor(getLabelColor_(coil_color_label));
				dp.setStrandColor(getLabelColor_(strand_color_label));
				dp.setTurnColor(getLabelColor_(turn_color_label));

				return;
			}

			if (RTTI::isKindOf<TemperatureFactorColorProcessor>(cp))
			{
				TemperatureFactorColorProcessor& dp = (*(TemperatureFactorColorProcessor*)&cp);
				dp.setMinColor(getLabelColor_(unassigned_tf_label));
				dp.getColors()[0] = (getLabelColor_(minimum_tf_label));
				dp.getColors()[1] = (getLabelColor_(maximum_tf_label));
				dp.setMaxColor(getLabelColor_(unassigned_tf_label));
				dp.setMaxValue(((float)max_tf_slider->value()) / 10.0);
				return;
			}

			if (RTTI::isKindOf<ForceColorProcessor>(cp))
			{
				ForceColorProcessor& dp = (*(ForceColorProcessor*)&cp);
				dp.getColors()[0] = (getLabelColor_(force_min_color_label));
				dp.getColors()[1] = (getLabelColor_(force_max_color_label));
				dp.setMaxValue(((float)force_max_value_slider->value()) / 10.0);
				dp.setMinValue(((float)force_min_value_slider->value()) / 10.0);
				return;
			}

			if (RTTI::isKindOf<ResidueTypeColorProcessor>(cp))
			{
				ResidueTypeColorProcessor& dp = (*(ResidueTypeColorProcessor*)&cp);
				dp.setBasicColor(getLabelColor_(basic_color_label));
				dp.setAcidicColor(getLabelColor_(acidic_color_label));
				dp.setAromaticColor(getLabelColor_(aromatic_color_label));
				dp.setPolarColor(getLabelColor_(polar_color_label));
				dp.setHydrophobicColor(getLabelColor_(hydrophobic_color_label));
				dp.setOtherColor(getLabelColor_(other_color_label));
				return;
			}

			if (RTTI::isKindOf<ChainColorProcessor>(cp))
			{
				((ChainColorProcessor*)&cp)->setColors(getColors(COLORING_CHAIN));
				return;
			}

			if (RTTI::isKindOf<MoleculeColorProcessor>(cp))
			{
				((MoleculeColorProcessor*)&cp)->setColors(getColors(COLORING_MOLECULE));
				return;
			}

		}
			
		void ColoringSettingsDialog::maxDistanceChanged()
		{
			String text = String(((float)max_distance_slider->value()) / 10.0);
			text = text.trimRight("0");
			if (text.hasSuffix(".")) text += "0";
			max_distance_value_label->setText(text.c_str());
		}
			
		void ColoringSettingsDialog::maxTFChanged()
		{
			String text = String(((float)max_tf_slider->value()) / 10.0);
			text = text.trimRight("0");
			if (text.hasSuffix(".")) text += "0";
			max_tf_label->setText(text.c_str());
		}

		void ColoringSettingsDialog::forceMaxValueChanged()
		{
			String text = String(((float)force_max_value_slider->value()) / 10.0);
			text = text.trimRight("0");
			if (text.hasSuffix(".")) text += "0";
			force_max_value_label->setText(text.c_str());
		}

		void ColoringSettingsDialog::forceMinValueChanged()
		{
			String text = String(((float)force_min_value_slider->value()) / 10.0);
			text = text.trimRight("0");
			if (text.hasSuffix(".")) text += "0";
			force_min_value_label->setText(text.c_str());
		}


		ColorProcessor* ColoringSettingsDialog::createColorProcessor(ColoringMethod method) const
			throw(Exception::InvalidOption)
		{
			ColorProcessor* color_processor = 0;

			switch(method)
			{
				case COLORING_ELEMENT:
					color_processor = new ElementColorProcessor;
					break;

				case COLORING_RESIDUE_NAME:
					color_processor = new ResidueNameColorProcessor;
					break;

				case COLORING_RESIDUE_INDEX:
					color_processor = new ResidueNumberColorProcessor;
					break;

				case COLORING_SECONDARY_STRUCTURE:
					color_processor = new SecondaryStructureColorProcessor;
					break;

				case COLORING_ATOM_CHARGE:
					color_processor = new AtomChargeColorProcessor;
					break;

				case COLORING_CUSTOM:
					color_processor = new CustomColorProcessor;
					break;

				case COLORING_DISTANCE:
					color_processor = new AtomDistanceColorProcessor;
					break;

				case COLORING_TEMPERATURE_FACTOR:
					color_processor = new TemperatureFactorColorProcessor;
					break;

				case COLORING_OCCUPANCY:
					color_processor = new OccupancyColorProcessor;
					break;

				case COLORING_FORCES:
					color_processor = new ForceColorProcessor;
					break;

				case COLORING_RESIDUE_TYPE:
					color_processor = new ResidueTypeColorProcessor;
					break;

				case COLORING_CHAIN:
					color_processor = new ChainColorProcessor;
					break;

				case COLORING_MOLECULE:
					color_processor = new MoleculeColorProcessor;
					break;

				default:
					throw(Exception::InvalidOption(__FILE__, __LINE__, method));
			}

			applySettingsTo(*color_processor);

			return color_processor;
		}

		void ColoringSettingsDialog::getSettings(const ColorProcessor& cp)
			throw()
		{

			if (RTTI::isKindOf<CustomColorProcessor>(cp))
			{
			} else

			if (RTTI::isKindOf<ElementColorProcessor>(cp))
			{
			} else
			
			if (RTTI::isKindOf<ResidueNameColorProcessor>(cp))
			{
			} else

			if (RTTI::isKindOf<ResidueNumberColorProcessor>(cp))
			{
				ResidueNumberColorProcessor& dp = (*(ResidueNumberColorProcessor*)&cp);
				setLabelColor_(first_residue_label, dp.getFirstColor());
				setLabelColor_(middle_residue_label, dp.getMiddleColor());
				setLabelColor_(last_residue_label, dp.getLastColor());
			} else

			if (RTTI::isKindOf<AtomChargeColorProcessor>(cp))
			{
				AtomChargeColorProcessor& dp = (*(AtomChargeColorProcessor*)&cp);
				setLabelColor_(negative_charge_label, dp.getColors()[0]);
				setLabelColor_(neutral_charge_label, dp.getColors()[1]);
				setLabelColor_(positive_charge_label, dp.getColors()[2]);
			} else

			if (RTTI::isKindOf<AtomDistanceColorProcessor>(cp))
			{
				AtomDistanceColorProcessor& dp = (*(AtomDistanceColorProcessor*)&cp);
				setLabelColor_(null_distance_label, dp.getNullDistanceColor());
				setLabelColor_(max_distance_label, dp.getMaxDistanceColor());
				max_distance_slider->setValue((Size)(dp.getDistance() * 10.0));
 				distance_show_selected->setChecked(dp.showSelected());
			} else

			if (RTTI::isKindOf<OccupancyColorProcessor>(cp))
			{
				OccupancyColorProcessor& dp = (*(OccupancyColorProcessor*)&cp);
				setLabelColor_(minimum_o_label, dp.getColors()[0]);
				setLabelColor_(maximum_o_label, dp.getColors()[1]);
			} else

			if (RTTI::isKindOf<SecondaryStructureColorProcessor>(cp))
			{
				SecondaryStructureColorProcessor& dp = (*(SecondaryStructureColorProcessor*)&cp);
				setLabelColor_(helix_color_label, dp.getHelixColor());
				setLabelColor_(coil_color_label, dp.getCoilColor());
				setLabelColor_(strand_color_label, dp.getStrandColor());
				setLabelColor_(turn_color_label, dp.getTurnColor());
			} else

			if (RTTI::isKindOf<TemperatureFactorColorProcessor>(cp))
			{
				TemperatureFactorColorProcessor& dp = (*(TemperatureFactorColorProcessor*)&cp);
				setLabelColor_(unassigned_tf_label, dp.getDefaultColor());
				setLabelColor_(minimum_tf_label, dp.getColors()[0]);
				setLabelColor_(maximum_tf_label, dp.getColors()[1]);
				max_tf_slider->setValue((Size)(dp.getMaxValue() * 10.0));
			} else

			if (RTTI::isKindOf<ForceColorProcessor>(cp))
			{
				ForceColorProcessor& dp = (*(ForceColorProcessor*)&cp);
				setLabelColor_(force_min_color_label, dp.getColors()[0]);
				setLabelColor_(force_max_color_label, dp.getColors()[1]);
				force_max_value_slider->setValue((Size)(dp.getMaxValue() * 10.0));
				force_min_value_slider->setValue((Size)(dp.getMinValue() * 10.0));
			} else

			if (RTTI::isKindOf<ResidueTypeColorProcessor>(cp))
			{
 				ResidueTypeColorProcessor& dp = (*(ResidueTypeColorProcessor*)&cp);
				setLabelColor_(acidic_color_label, dp.getAcidicColor());
				setLabelColor_(aromatic_color_label, dp.getAromaticColor());
				setLabelColor_(basic_color_label, dp.getBasicColor());
				setLabelColor_(hydrophobic_color_label, dp.getHydrophobicColor());
				setLabelColor_(other_color_label, dp.getOtherColor());
				setLabelColor_(polar_color_label, dp.getPolarColor());
			} else

			if (RTTI::isKindOf<ChainColorProcessor>(cp))
			{
 				ChainColorProcessor& dp = (*(ChainColorProcessor*)&cp);
				vector<String> 		names;
				vector<ColorRGBA> colors;

				for (Position p = 0; p < dp.getColors().size(); p++)
				{
					colors.push_back(dp.getColors()[p]);
					names.push_back(p);
				}

				chain_table_->setNamesTitle("Chain");
				chain_table_->setContent(names, colors);
			} else

			if (RTTI::isKindOf<MoleculeColorProcessor>(cp))
			{
 				MoleculeColorProcessor& dp = (*(MoleculeColorProcessor*)&cp);
				vector<String> 		names;
				vector<ColorRGBA> colors;

				for (Position p = 0; p < dp.getColors().size(); p++)
				{
					colors.push_back(dp.getColors()[p]);
					names.push_back(p);
				}

				molecule_table_->setNamesTitle("Molecule");
				molecule_table_->setContent(names, colors);
			}

		}

		QWidget* ColoringSettingsDialog::getEntryFor(ColoringMethod method)
			throw()
		{
			switch (method)
			{
				case COLORING_ELEMENT: 							return widget_stack->widget(0);
				case COLORING_RESIDUE_NAME: 				return widget_stack->widget(2);
				case COLORING_RESIDUE_INDEX: 				return widget_stack->widget(1);
				case COLORING_SECONDARY_STRUCTURE: 	return widget_stack->widget(7);
				case COLORING_ATOM_CHARGE: 					return widget_stack->widget(3);
				case COLORING_DISTANCE: 						return widget_stack->widget(4);
				case COLORING_TEMPERATURE_FACTOR: 	return widget_stack->widget(5);
				case COLORING_OCCUPANCY: 						return widget_stack->widget(6);
				case COLORING_FORCES: 							return widget_stack->widget(8);
				case COLORING_RESIDUE_TYPE: 				return widget_stack->widget(9);
				case COLORING_CHAIN: 								return widget_stack->widget(10);
				default: break;
			}

			return 0;
		}

  } // namespace VIEW
} // namespace BALL
