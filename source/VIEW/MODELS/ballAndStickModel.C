// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// $Id: ballAndStickModel.C,v 1.22.2.1 2005/10/12 16:14:45 amoll Exp $
//

#include <BALL/VIEW/MODELS/ballAndStickModel.h>
#include <BALL/KERNEL/atom.h>
#include <BALL/KERNEL/bond.h>
#include <BALL/KERNEL/forEach.h>
#include <BALL/VIEW/PRIMITIVES/tube.h>
#include <BALL/VIEW/PRIMITIVES/sphere.h>
#include <BALL/VIEW/PRIMITIVES/disc.h>
#include <BALL/VIEW/PRIMITIVES/twoColoredTube.h>
#include <BALL/VIEW/KERNEL/common.h>

using namespace std;

namespace BALL
{
	namespace VIEW
	{

		AddBallAndStickModel::AddBallAndStickModel()
			throw()
			: AtomBondModelBaseProcessor(),
				ball_radius_((float)0.4),
				stick_radius_((float)0.2),
				ball_and_stick_(true),
				dashed_bonds_(false)
		{
		}

		AddBallAndStickModel::AddBallAndStickModel(const AddBallAndStickModel &add_ball_and_stick)
			throw()
			: AtomBondModelBaseProcessor(add_ball_and_stick),
				ball_radius_(add_ball_and_stick.ball_radius_),
				stick_radius_(add_ball_and_stick.stick_radius_),
				ball_and_stick_(add_ball_and_stick.ball_and_stick_),
				dashed_bonds_(add_ball_and_stick.dashed_bonds_)
		{
		}

		AddBallAndStickModel::~AddBallAndStickModel()
			throw()
		{
			#ifdef BALL_VIEW_DEBUG
				Log.info() << "Destructing object " << this << " of class AddBallAndStickModel" << std::endl;
			#endif 
		}

		void AddBallAndStickModel::clear()
			throw()
		{
			AtomBondModelBaseProcessor::clear();

			ball_radius_ = (float)0.4;
			stick_radius_ = (float)0.2;
			ball_and_stick_ = true;
			dashed_bonds_ = false;
		}

		void AddBallAndStickModel::set(const AddBallAndStickModel &add_ball_and_stick)
			throw()
		{
			AtomBondModelBaseProcessor::set(add_ball_and_stick);

			ball_radius_ = add_ball_and_stick.ball_radius_;
			stick_radius_ = add_ball_and_stick.stick_radius_;
			ball_and_stick_ = add_ball_and_stick.ball_and_stick_;
			dashed_bonds_ = add_ball_and_stick.dashed_bonds_;
		}

		const AddBallAndStickModel &AddBallAndStickModel::operator = (const AddBallAndStickModel &processor)
			throw()
		{
			set(processor);
			return *this;
		}

		void AddBallAndStickModel::setBallRadius(const float radius)
			throw(Exception::OutOfRange)
		{
			// a radius never can be lower or equal 0
			if (radius <= (float)0)
			{
				throw Exception::OutOfRange(__FILE__, __LINE__);
			}
			
			ball_radius_ = radius;
		}

		void AddBallAndStickModel::setStickRadius(const float radius)
			throw(Exception::OutOfRange)
		{
			// a radius never can be lower or equal 0
			if (radius <= (float)0)
			{
				throw Exception::OutOfRange(__FILE__, __LINE__);
			}

			stick_radius_ = radius;
		}

		Processor::Result AddBallAndStickModel::operator() (Composite& composite)
		{
			if (dashed_bonds_)
			{
				AtomBondModelBaseProcessor::operator() (composite);
			}

			Atom *atom = RTTI::castTo<Atom>(composite);

			if (atom == 0) return Processor::CONTINUE;

			Sphere* sphere_ptr = new Sphere;

			if (sphere_ptr == 0) throw Exception::OutOfMemory (__FILE__, __LINE__, sizeof(Sphere));

			sphere_ptr->setComposite(atom);

			if (ball_and_stick_)
			{
				sphere_ptr->setRadius(ball_radius_);
			}
			else
			{
				sphere_ptr->setRadius(stick_radius_);
			}

			sphere_ptr->setPositionAddress(atom->getPosition());
			
			// append sphere in Atom
			geometric_objects_.push_back(sphere_ptr);

			// collect used atoms
			insertAtom_(atom);

			return Processor::CONTINUE;
		}

		void AddBallAndStickModel::dump(std::ostream& s, Size depth) const
			throw()
		{
			BALL_DUMP_STREAM_PREFIX(s);
			
			BALL_DUMP_DEPTH(s, depth);
			BALL_DUMP_HEADER(s, this, this);

			AtomBondModelBaseProcessor::dump(s, depth + 1);

			BALL_DUMP_DEPTH(s, depth);
			s << "ball radius: " << ball_radius_ << endl;

			BALL_DUMP_DEPTH(s, depth);
			s << "stick radius: " << stick_radius_ << endl;

			BALL_DUMP_DEPTH(s, depth);
			s << "b&s model: " << ball_and_stick_ << endl;

			BALL_DUMP_DEPTH(s, depth);
			s << "s model: " << ball_and_stick_ << endl;

			BALL_DUMP_STREAM_SUFFIX(s);
		}

		void AddBallAndStickModel::visualiseBond_(const Bond& bond)
			throw()
		{
			// no visualisation for hydrogen bonds
			if (bond.getType() == Bond::TYPE__HYDROGEN) return;

			if (!dashed_bonds_ 													||
					bond.getOrder() < Bond::ORDER__DOUBLE 	||
					bond.getOrder() > Bond::ORDER__AROMATIC )
			{
				// generate two colored tube
				TwoColoredTube *tube = new TwoColoredTube;
				tube->setRadius(stick_radius_);
				tube->setVertex1Address(bond.getFirstAtom()->getPosition());
				tube->setVertex2Address(bond.getSecondAtom()->getPosition());
				tube->setComposite(&bond);
				geometric_objects_.push_back(tube);
				return;
			}

			try
			{
				Vector3 dir = bond.getSecondAtom()->getPosition() - bond.getFirstAtom()->getPosition();
				Vector3 normal = VIEW::getNormal(dir);

				if (bond.getOrder() == Bond::ORDER__AROMATIC)
				{
					// Bonds in aromatic rings will be drawn later
					if (ring_atoms_.has(bond.getFirstAtom()) &&
							ring_atoms_.has(bond.getSecondAtom()))
					{
						return;
					}

					renderDashedBond_(*bond.getFirstAtom(), *bond.getSecondAtom(), normal, normal);
					return;
				}
	
				normal *= stick_radius_ / (float) 1.5;
				float radius = stick_radius_ / (float) 2.4;

				TwoColoredTube* tube = new TwoColoredTube;
				tube->setRadius(radius);
				tube->setComposite(&bond);
				geometric_objects_.push_back(tube);
					
				TwoColoredTube* tube2 = new TwoColoredTube(*tube);
				geometric_objects_.push_back(tube2);

				if (bond.getOrder() == Bond::ORDER__DOUBLE)
				{
					tube->setVertex1(bond.getFirstAtom()->getPosition() - normal);
					tube->setVertex2(bond.getSecondAtom()->getPosition() - normal);

					tube2->setVertex1(bond.getFirstAtom()->getPosition() + normal);
					tube2->setVertex2(bond.getSecondAtom()->getPosition() + normal);
				}
				else
				{
					normal *= radius;

					Vector3 normal2 = dir % normal;
					normal2.normalize();
					normal2 *= radius;
					
					tube->setVertex1(bond.getFirstAtom()->getPosition() - normal - normal2);
					tube->setVertex2(bond.getSecondAtom()->getPosition() - normal - normal2);
					
					tube2->setVertex1(bond.getFirstAtom()->getPosition() + normal - normal2);
					tube2->setVertex2(bond.getSecondAtom()->getPosition() + normal - normal2);

					TwoColoredTube* tube3 = new TwoColoredTube(*tube);
					tube3->setVertex1(bond.getFirstAtom()->getPosition() + normal2);
					tube3->setVertex2(bond.getSecondAtom()->getPosition() + normal2);
					geometric_objects_.push_back(tube3);
				}
			}
			catch(...)
			{
				Log.error() << "We have a bond with a length of 0. Can't draw this!" << std::endl;
			}
		}

		void AddBallAndStickModel::visualiseRings_()
			throw()
		{
			if (!dashed_bonds_) return;

			vector<vector<Atom*> >::iterator it = rings_.begin();
			for(; it != rings_.end(); it++)
			{
				vector<Atom*>& ring = *it;
				if (ring.size() != 5 && ring.size() != 6) continue;

				vector<Atom*>::iterator ait = ring.begin();

				Vector3 center;
				for (; ait != ring.end(); ait++)
				{
					center += (**ait).getPosition();
				}
				center /= ring.size();

				for (ait = ring.begin(); ait != ring.end(); ait++)
				{
					AtomBondIterator bi;
					BALL_FOREACH_ATOM_BOND(**ait, bi)
					{
						Atom* partner = bi->getPartner(**ait);
						if (partner < *ait) continue;
						
						Vector3 n1(center - (**ait).getPosition());
						Vector3 n2(center - partner->getPosition());

						renderDashedBond_(**ait, *partner, n1, n2);
					}
				}
			}
		}

		void AddBallAndStickModel::renderDashedBond_(const Atom& a1, const Atom& a2, Vector3 n1, Vector3 n2)
			throw(Exception::DivisionByZero)
		{
			const Bond& bond = *a1.getBond(a2);

			if (bond.getOrder() != Bond::ORDER__AROMATIC) return;

			n1.normalize();
			n2.normalize();
			n1 *= stick_radius_ / (float) 1.5;
			n2 *= stick_radius_ / (float) 1.5;

			// render one tube with full length			
			TwoColoredTube *tube = new TwoColoredTube;
			tube->setRadius(stick_radius_ / (float) 2.4);
			tube->setVertex1(a1.getPosition() - n1);
			tube->setVertex2(a2.getPosition() - n2);
			tube->setComposite(&bond);
			geometric_objects_.push_back(tube);

			// render dashed tubes
			Vector3 v = a2.getPosition() + n2 - (a1.getPosition() + n1);
			Vector3 last = a1.getPosition() + n1 + v / (float) 4.5;
			for (Position p = 0; p < 3; p++)
			{
				if (p == 1)
				{
					Vector3 middle((a1.getPosition() - a2.getPosition()) / 2 + a2.getPosition());
					TwoColoredTube *tube = new TwoColoredTube;
					tube->setRadius(stick_radius_ / (float) 2.4);
					tube->setComposite(&bond);
					tube->setVertex1(middle - (v / 8) + n1);
					tube->setVertex2(middle + (v / 8) + n2);
					geometric_objects_.push_back(tube);

					Disc* disc = new Disc(Circle3(middle - (v / 8) + n1, -v, stick_radius_ / (float) 2.4));
					disc->setComposite(&a1);
					geometric_objects_.push_back(disc);

					disc = new Disc(Circle3(middle + (v / 8) + n2, v, stick_radius_ / (float) 2.4));
					disc->setComposite(&a2);
					geometric_objects_.push_back(disc);

					last += (v /4);

					continue;
				}

				last += (v /4);
			}
		}

#		ifdef BALL_NO_INLINE_FUNCTIONS
#			include <BALL/VIEW/MODELS/ballAndStickModel.iC>
#		endif

	} // namespace VIEW
} // namespace BALL
