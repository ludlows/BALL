// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// $Id: colorProcessor.C,v 1.34.4.12 2005/11/24 14:17:06 amoll Exp $
//

#include <BALL/VIEW/MODELS/colorProcessor.h>

#include <BALL/COMMON/limits.h>
#include <BALL/VIEW/DATATYPE/colorExtension2.h>
#include <BALL/VIEW/PRIMITIVES/mesh.h>
#include <BALL/KERNEL/bond.h>
#include <BALL/KERNEL/forEach.h>
#include <BALL/KERNEL/PTE.h>
#include <BALL/DATATYPE/list.h>
#include <BALL/STRUCTURE/geometricProperties.h>
#include <BALL/SYSTEM/sysinfo.h>

using namespace std;

namespace BALL
{
	namespace VIEW
	{

		ColorProcessor::ColorProcessor()
			throw()
			:	UnaryProcessor<GeometricObject*>(),
				update_always_needed_(false)
		{
			clear();
		}

		ColorProcessor::ColorProcessor(const ColorProcessor& cp)
			throw()
			:	UnaryProcessor<GeometricObject*>(cp),
				update_always_needed_(cp.update_always_needed_),
				default_color_(cp.default_color_),
				selection_color_(cp.selection_color_),
				transparency_(cp.transparency_),
				composites_(cp.composites_),
				atom_grid_(cp.atom_grid_),
				model_type_(cp.model_type_),
				last_composite_of_grid_(cp.last_composite_of_grid_),
				additional_grid_distance_(cp.additional_grid_distance_)
		{
		}

		ColorProcessor::~ColorProcessor()
		{
			#ifdef BALL_VIEW_DEBUG
				Log.error() << "Destructing object " << (void *)this 
										<< " of class " << RTTI::getName<ColorProcessor>() << std::endl;
			#endif 
		}

		void ColorProcessor::clear()
			throw()
		{
			default_color_.set("FF0000FF");
			transparency_ = 0;
			composites_ = 0;
			clearAtomGrid();
			selection_color_ = BALL_SELECTED_COLOR;
			last_composite_of_grid_ = 0;
			additional_grid_distance_ = 5.0;
		}

		void ColorProcessor::set(const ColorProcessor& cp)
			throw()
		{
			update_always_needed_ = cp.update_always_needed_;
			default_color_ = cp.default_color_;
			selection_color_ = cp.selection_color_;
			transparency_ = cp.transparency_;
			composites_ = cp.composites_;
			atom_grid_ = cp.atom_grid_;
			model_type_ = cp.model_type_;
			last_composite_of_grid_ = cp.last_composite_of_grid_;
			additional_grid_distance_ = cp.additional_grid_distance_;
		}


		const ColorProcessor& ColorProcessor::operator = (const ColorProcessor& cp)
			throw()
		{
			set(cp);
			return *this;
		}


		void ColorProcessor::dump(ostream& s, Size depth) const
			throw()
		{
			BALL_DUMP_STREAM_PREFIX(s);
			
			BALL_DUMP_DEPTH(s, depth);
			BALL_DUMP_HEADER(s, this, this);

			BALL_DUMP_DEPTH(s, depth);
			s << "default_color: " << default_color_ << endl;
					
			BALL_DUMP_STREAM_SUFFIX(s);
		}

		Processor::Result ColorProcessor::operator() (GeometricObject*& object)
		{
			const Composite* composite = object->getComposite();

			/// ------------- meshes ----------------------------
			Mesh* const mesh = dynamic_cast<Mesh*>(object);
			if (mesh != 0)
			{
				mesh->colors.clear();
				if (composite == &composite_to_be_ignored_for_colorprocessors_ ||
						composites_ == 0)
				{
					mesh->colors.push_back(default_color_);
					return Processor::CONTINUE;
				}

				if (composite == 0 ||
						composite != last_composite_of_grid_)
				{
					createAtomGrid(composite);
				}

				colorMeshFromGrid_(*mesh);
				return Processor::CONTINUE;
			}

			ColorExtension2* const two_colored = dynamic_cast<ColorExtension2*>(object);

			/// ------------- custom objects ----------------------------
			if (composite == 0 ||
					composite == &composite_to_be_ignored_for_colorprocessors_)
			{
				object->setColor(default_color_); 
				if (two_colored != 0)
				{
					two_colored->setColor2(default_color_);
				}
				return Processor::CONTINUE;
			}

			/// ------------- single colored objects ----------------------------
			if (two_colored == 0)
			{
				if (composite->isSelected())
				{
					object->setColor(selection_color_);
				}
				else
				{
					getColor(*composite, object->getColor()); 
				}
				return Processor::CONTINUE;
			}

			/// ------------- two colored objects ----------------------------
			const Bond* const bond = dynamic_cast<const Bond*>(composite);


			/// ------------- bonds  ----------------------------
			if (bond != 0)
			{
				const Atom* atom = bond->getFirstAtom();
				if (!atom->isSelected())
				{
					getColor(*atom, object->getColor());
				}
				else
				{
					object->setColor(selection_color_);
				}

				const Atom* atom2 = bond->getSecondAtom();
				if (!atom2->isSelected())
				{
					getColor(*atom2, two_colored->getColor2());
				}
				else
				{
					two_colored->setColor2(selection_color_);
				}
			}
			/// ------------- non bonds  ----------------------------
			else
			{
				if (composite->isSelected())
				{
					object->setColor(selection_color_);
					two_colored->setColor2(selection_color_);
				}
				else
				{
 					getColor(*composite, object->getColor());
 					two_colored->setColor2(object->getColor());
				}
			}

			return Processor::CONTINUE;
		}

		void ColorProcessor::createAtomGrid(const Composite* from_mesh)
			throw()
		{
			atom_grid_.clear();
			if (composites_ == 0 && from_mesh == 0)
			{
				return;
			}
			
			List<const Atom*> atoms;

			if (from_mesh == 0)
			{
				List<const Composite*>::const_iterator it = composites_->begin();
				for(; it != composites_->end(); it++)
				{
					if (RTTI::isKindOf<AtomContainer>(**it))
					{
						AtomConstIterator ait;
						const AtomContainer* const acont = dynamic_cast<const AtomContainer*>(*it);
						BALL_FOREACH_ATOM(*acont, ait)
						{
							atoms.push_back(&*ait);
						}
					}
					else if (RTTI::isKindOf<Atom>(**it))
					{
						const Atom* atom = dynamic_cast<const Atom*> (*it);
						atoms.push_back(atom);
					}
				}
			}
			else 
			{
				// composite from mesh
				if (RTTI::isKindOf<AtomContainer>(*from_mesh))
				{
					AtomConstIterator ait;
					const AtomContainer* const acont = dynamic_cast<const AtomContainer*>(from_mesh);
					BALL_FOREACH_ATOM(*acont, ait)
					{
						atoms.push_back(&*ait);
					}
				}
				else if (RTTI::isKindOf<Atom>(*from_mesh))
				{
					const Atom* atom = dynamic_cast<const Atom*> (from_mesh);
					atoms.push_back(atom);
				}
			}


			BoundingBoxProcessor boxp;
			boxp.start();
			List<const Atom*>::Iterator lit = atoms.begin();
			for(;lit != atoms.end(); lit++)
			{
				boxp.operator() (*(Atom*)*lit);
			}
			boxp.finish();

			const Vector3 diagonal = boxp.getUpper() - boxp.getLower();
			
			// abort for strange molecules with huge dimension
			if (diagonal.getSquareLength() > 10000000.0) return;

			// grid spacing, tradeoff between speed and memory consumption
			float grid_spacing = 4.0;
			if (diagonal.getSquareLength() < 5000)
			{
				grid_spacing = 3.5;
			} 
			else 
			{
				float memory = SysInfo::getAvailableMemory();
				// if we can not calculate available memory, use around 60 MB for the grid
				if (memory == -1) memory = 10000000;
				memory *= 0.6;
				float min_spacing = HashGrid3<const Atom*>::calculateMinSpacing((LongIndex)memory, diagonal + 
																																						Vector3(2 * (additional_grid_distance_ + 15.0)));
				if (min_spacing > grid_spacing) grid_spacing = min_spacing;
			}
			
			atom_grid_ = AtomGrid(boxp.getLower() - Vector3(additional_grid_distance_ + 15.0),
														diagonal + Vector3(2 * additional_grid_distance_ + 15.0),
														grid_spacing); 
		 
			for (lit = atoms.begin(); lit != atoms.end(); lit++)
			{
				atom_grid_.insert((*lit)->getPosition(), *lit);
			}

			last_composite_of_grid_ = from_mesh;
		}

		void ColorProcessor::colorMeshFromGrid_(Mesh& mesh)
			throw()
		{
			if (atom_grid_.isEmpty()) return;
			
			mesh.colors.resize(mesh.vertex.size());
			
			for (Position p = 0; p < mesh.vertex.size(); p++)
			{
				// make sure we found an atom
				const Atom* atom = getClosestItem(mesh.vertex[p]);

				if (atom == 0)
				{
					mesh.colors[p] = default_color_;
				}
				else
				{
					if (atom->isSelected())
					{
						mesh.colors[p] = selection_color_;
					}
					else
					{
						getColor(*atom, mesh.colors[p]);
					}
				}
			}
		}

		void ColorProcessor::setComposites(const List<const Composite*>* composites)
			throw() 
		{ 
			composites_ = composites;
			clearAtomGrid();
		}

		void ColorProcessor::clearAtomGrid()
			throw()
		{
			atom_grid_.clear();
		}

		void ColorProcessor::setDefaultColor(const ColorRGBA& color)
			throw()
		{
			default_color_ = color;
			default_color_.setAlpha(255 - transparency_);
		}

		void ColorProcessor::setTransparency(Size value)
			throw() 
		{ 
			transparency_ = value;
			default_color_.setAlpha(255 - value);
			selection_color_.setAlpha(255 - value);
		}

		void ColorProcessor::getColor(const Composite& composite, ColorRGBA& color_to_be_set)
		{
			if (composite.isSelected())
			{
				color_to_be_set = selection_color_;
			}
			else
			{
				color_to_be_set = default_color_;
			}
		}


		const Atom* ColorProcessor::getClosestItem(const Vector3& point) const
			throw()
		{
			const HashGridBox3<const Atom*>* box = atom_grid_.getBox(point);
			if (!box)
			{
				return 0;
			}

			Position x, y, z;
			atom_grid_.getIndices(*box, x, y, z);

			const Atom* const* item = 0;
			float distance = Limits<float>::max();
			List<HashGridBox3<const Atom*>* > box_list;
			Size dist = 1;
			// iterator over neighbour boxes
			for (Index xi = -(Index)dist; xi <= (Index)dist; xi++)
			{
				for (Index yi = -(Index)dist; yi <= (Index)dist; yi++)
				{
					for (Index zi = -(Index)dist; zi <= (Index)dist; zi++)
					{
						// iterate over all data items
						const HashGridBox3<const Atom*>* box_ptr = atom_grid_.getBox(x+xi, y+yi, z+zi);	
						if (box_ptr != 0 && !box_ptr->isEmpty())
						{
							HashGridBox3<const Atom*>::ConstDataIterator hit = box_ptr->beginData();
							for (;hit != box_ptr->endData(); hit++)
							{
								// this is not 
								float radius = (*hit)->getElement().getVanDerWaalsRadius();
								if (model_type_ == MODEL_SA_SURFACE)
								{
									radius += 2;
								}
								if (radius <= 0.0) radius = 1;
								// avoid calculation of the square roots
								float new_dist = ((*hit)->getPosition() - point).getSquareLength() - radius * radius;
								if (new_dist < distance)
								{
									item = &*hit;
									distance = new_dist;
								}
							}
						}
					}
				}
			}

			if (item == 0) return 0;
			
			return *item;
		}

		bool ColorProcessor::start()
			throw()
		{
			selection_color_ = BALL_SELECTED_COLOR;
			selection_color_.setAlpha(255 - transparency_);
			return (getComposites() != 0);
		}

		//////////////////////////////////////////////////////////////////////
		InterpolateColorProcessor::InterpolateColorProcessor()
			: ColorProcessor(),
				min_color_(ColorRGBA(1.0,1.0,1.0)),
				max_color_(ColorRGBA(1.0,0.0,0)),
				mode_(NO_OUTSIDE_COLORS),
				max_value_(1),
				min_value_(0)
		{
			default_color_ = ColorRGBA(1.0,1.0,1.0);
		}


		InterpolateColorProcessor::InterpolateColorProcessor(const InterpolateColorProcessor& pro)
			: ColorProcessor(pro)
		{
		}

		
		bool InterpolateColorProcessor::start()
			throw()
		{
			if (!ColorProcessor::start()) return false;

			if (colors_.size() < 2 ||
					max_value_ <= min_value_)
			{
				return false;
			}

			x_ = (max_value_ - min_value_) / (colors_.size() - 1);

			min_color_.setAlpha(255 - transparency_);
			max_color_.setAlpha(255 - transparency_);
			default_color_.setAlpha(255 - transparency_);

			for (Position p = 0; p < colors_.size(); p++)
			{
				colors_[p].setAlpha(255 - transparency_);
			}

			if (mode_ == NO_OUTSIDE_COLORS)
			{
				min_color_ = colors_[0];
				max_color_ = colors_[colors_.size() - 1];
			}
			else if (mode_ == DEFAULT_COLOR_FOR_OUTSIDE_COLORS)
			{
				min_color_ = default_color_;
				max_color_ = default_color_;
			}

			return true;
		}

		void InterpolateColorProcessor::interpolateColor(float value, ColorRGBA& color_to_be_set)
			throw()
		{
			if (value < min_value_)
			{
				color_to_be_set.set(min_color_);
				return;
			}
			
			if (value > max_value_) 
			{
				color_to_be_set.set(max_color_);
				return;
			}

			const Position z1 = (Position)floor((value - min_value_)/ x_);
 			const Position z2 = (Position)ceil((value - min_value_)/ x_);

			const float& red1   = colors_[z1].getRed();
			const float& green1 = colors_[z1].getGreen();
			const float& blue1  = colors_[z1].getBlue();

			const float red2   = (float) colors_[z2].getRed()   - (float) colors_[z1].getRed();
			const float green2 = (float) colors_[z2].getGreen() - (float) colors_[z1].getGreen();
			const float blue2  = (float) colors_[z2].getBlue()  - (float) colors_[z1].getBlue();

			const float dz1 = (value - min_value_ - ((float)z1) * x_) / x_;

			color_to_be_set.set(red1   + dz1 * red2,
													green1 + dz1 * green2,
													blue1  + dz1 * blue2,
													255 - transparency_);
		}

		void InterpolateColorProcessor::setMinColor(const ColorRGBA& color)
			throw()
		{
			min_color_ = color;
		}

		void InterpolateColorProcessor::setMaxColor(const ColorRGBA& color)
			throw()
		{
			max_color_ = color;
		}

		const ColorRGBA& InterpolateColorProcessor::getMinColor() const
			throw()
		{
			return min_color_;
		}

		const ColorRGBA& InterpolateColorProcessor::getMaxColor() const
			throw()
		{
			return max_color_;
		}

	} // namespace VIEW
 } // namespace BALL
