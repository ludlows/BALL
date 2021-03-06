// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// $Id: buildBondsProcessor.C,v 1.7.2.4 2005/05/01 11:37:22 oliver Exp $
//

#include <BALL/STRUCTURE/buildBondsProcessor.h>
#include <BALL/KERNEL/PTE.h>
#include <BALL/KERNEL/forEach.h>
#include <BALL/KERNEL/atom.h>
#include <BALL/KERNEL/bond.h>
#include <BALL/DATATYPE/hashGrid.h>
#include <BALL/COMMON/limits.h>
#include <BALL/SYSTEM/path.h>
#include <BALL/FORMAT/resourceFile.h>
#include <BALL/STRUCTURE/geometricProperties.h>
#include <BALL/QSAR/ringPerceptionProcessor.h>

using namespace std;

namespace BALL 
{
	const char* BuildBondsProcessor::Option::BONDLENGTHS_FILENAME = "bondlengths_filename";
	const char* BuildBondsProcessor::Option::DELETE_EXISTING_BONDS = "delete_existing_bonds";
	const char* BuildBondsProcessor::Option::REESTIMATE_BONDORDERS_RINGS = "reestimate_bondorders_rings";
	const char* BuildBondsProcessor::Option::DELETE_OVERESTIMATED_BONDS = "delete_overestimated_bonds";
	const char* BuildBondsProcessor::Default::BONDLENGTHS_FILENAME = "bond_lengths/bond_lengths.db";
	const bool  BuildBondsProcessor::Default::DELETE_EXISTING_BONDS = false;
	const bool  BuildBondsProcessor::Default::REESTIMATE_BONDORDERS_RINGS = false;
	const bool  BuildBondsProcessor::Default::DELETE_OVERESTIMATED_BONDS = false;
	
	BuildBondsProcessor::BuildBondsProcessor()
		: UnaryProcessor<AtomContainer>(),
			options(),
			num_bonds_(0),
			max_length_(0.0f)
	{
		readBondLengthsFromFile_();
	}

	BuildBondsProcessor::BuildBondsProcessor(const BuildBondsProcessor& bbp)
		:	UnaryProcessor<AtomContainer>(bbp),
			options(),
			num_bonds_(0),
			max_length_(bbp.max_length_)
	{
		bond_lengths_ = bbp.bond_lengths_;
		max_bond_lengths_ = bbp.max_bond_lengths_;
		min_bond_lengths_ = bbp.min_bond_lengths_;
	}

	BuildBondsProcessor::BuildBondsProcessor(const String& file_name)	throw(Exception::FileNotFound)
		:	UnaryProcessor<AtomContainer>(),
			options(),
			num_bonds_(0),
			max_length_(0.0f)
	{
		readBondLengthsFromFile_(file_name);
	}

	BuildBondsProcessor::~BuildBondsProcessor()
	{
	}

	BuildBondsProcessor& BuildBondsProcessor::operator = (const BuildBondsProcessor& bbp)
	{
		max_length_ = bbp.max_length_;
		bond_lengths_ = bbp.bond_lengths_;
		max_bond_lengths_ = bbp.max_bond_lengths_;
		min_bond_lengths_ = bbp.min_bond_lengths_;
		return *this;
	}

	bool BuildBondsProcessor::start()
	{
		num_bonds_ = 0;
		return true;
	}

	Size BuildBondsProcessor::getNumberOfBondsBuilt()
	{
		return num_bonds_;
	}


	Processor::Result BuildBondsProcessor::operator () (AtomContainer& ac)
	{
		if (options.getBool(BuildBondsProcessor::Option::DELETE_EXISTING_BONDS))
		{
			AtomIterator ait;
			BALL_FOREACH_ATOM(ac, ait)
			{
				ait->destroyBonds();
			}
		}
		
		num_bonds_ += buildBondsHashGrid3_(ac);

		estimateBondOrders_(ac);

		if (options.getBool(BuildBondsProcessor::Option::REESTIMATE_BONDORDERS_RINGS))
		{
			reestimateBondOrdersRings_(ac);
		}
	
		if (options.getBool(BuildBondsProcessor::Option::DELETE_OVERESTIMATED_BONDS))
		{
			deleteOverestimatedBonds_(ac);
		}
		
		return Processor::BREAK;
	}

	Size BuildBondsProcessor::buildBondsHashGrid3_(AtomContainer& ac)
	{
		Size num_bonds(0);

		// get the bounding box of ac
		BoundingBoxProcessor bbox;
		ac.apply(bbox);		
		Vector3 size = bbox.getUpper() - bbox.getLower();

		// build HashGrid
		HashGrid3<Atom*> grid(bbox.getLower() - Vector3(1.0), size + Vector3(2.0), max_length_ + 0.01);
	
		// fill in the atom pointers into the grid
		for (AtomIterator a_it(ac.beginAtom()); +a_it; ++a_it)
		{
			grid.insert(a_it->getPosition(), &*a_it);
		}
	

		// iterate over all boxes
		HashGrid3<Atom*>::BoxIterator box_it = grid.beginBox();
		for (; +box_it; ++box_it) 
		{
			Position x, y, z;
			grid.getIndices(*box_it, x, y, z);

			// iterator over neighbour boxes
			for (Index xi = -1; xi <= 1; xi++)
			{
				const Position nx = x + xi;
				
				for (Index yi = -1; yi <= 1; yi++)
				{
					const Position ny = y + yi;

					for (Index zi = -1; zi <= 1; zi++)
					{
						HashGridBox3<Atom*>* const bbox = grid.getBox(nx, ny, z+zi);
						// smaller operator also checks for 0 !
						if (bbox < &*box_it || bbox->isEmpty()) 
						{
							continue;
						}

						if (bbox == &*box_it)
						{
							HashGridBox3<Atom*>::ConstDataIterator ait1 = box_it->beginData();
							for (;+ait1; ait1++)
							{
								const Vector3& atom_pos1 = (*ait1)->getPosition();
								Atom& atom1 = **ait1;

								HashGridBox3<Atom*>::ConstDataIterator ait2 = bbox->beginData();
								for (;+ait2; ait2++)
								{
									if (*ait1 <= *ait2) continue;
									// test the distance criterion for the specific element pair
									const float dist = atom_pos1.getSquareDistance((*ait2)->getPosition());
									float max_dist, min_dist;
									const Size& an1 = atom1.getElement().getAtomicNumber();
									const Size& an2 = (*ait2)->getElement().getAtomicNumber();

									if (getMaxBondLength_(max_dist, an1, an2) &&
											getMinBondLength_(min_dist, an1, an2) &&
											dist <= max_dist && 
											dist >= min_dist &&
											!atom1.isBoundTo(**ait2))
									{
										atom1.createBond(**ait2);
										num_bonds++;
									}
								}
							}
							continue;
						}

						// iterate over all atoms of current box
						HashGridBox3<Atom*>::ConstDataIterator ait1 = box_it->beginData();
						for (;+ait1; ait1++)
						{
							const Vector3& atom_pos1 = (*ait1)->getPosition();
							Atom& atom1 = **ait1;

							HashGridBox3<Atom*>::ConstDataIterator ait2 = bbox->beginData();
							for (;+ait2; ait2++)
							{
								// test the distance criterion for the specific element pair
								const float dist = atom_pos1.getSquareDistance((*ait2)->getPosition());
								float max_dist, min_dist;
								const Size& an1 = atom1.getElement().getAtomicNumber();
								const Size& an2 = (*ait2)->getElement().getAtomicNumber();

								if (getMaxBondLength_(max_dist, an1, an2) &&
										getMinBondLength_(min_dist, an1, an2) &&
										dist <= max_dist && 
										dist >= min_dist &&
										!atom1.isBoundTo(**ait2))
								{
									Bond* const b = atom1.createBond(**ait2);
									b->setOrder(Bond::ORDER__UNKNOWN);
									num_bonds++;
								}
							} // data iterator ait2
						} // data iterator ait1
					} // zi
				} // yi
			}  // xi
		} // for all boxes

		return num_bonds;
	}
	
	void BuildBondsProcessor::estimateBondOrders_(AtomContainer& ac)
	{
		// iterate over all bonds
		AtomIterator ait;
		Atom::BondIterator bit;
		BALL_FOREACH_BOND(ac, ait, bit)
		{
			// set best bond order found
			if (bit->getOrder() == Bond::ORDER__UNKNOWN)
			{
				float length = bit->getLength();
				Size an1(bit->getFirstAtom()->getElement().getAtomicNumber());
				Size an2(bit->getSecondAtom()->getElement().getAtomicNumber());
				
				bit->setOrder(getNearestBondOrder_(length, an1, an2));
			}
		}
	}

	void BuildBondsProcessor::reestimateBondOrdersRings_(AtomContainer& ac)
	{
		RingPerceptionProcessor rpp;
		vector<vector<Atom*> > sssr;

		// abort if no ring calculated
		if (rpp.calculateSSSR(sssr, ac) == 0) return;
	
		vector<vector<Atom*> >::iterator it = sssr.begin();
		for (; it!=sssr.end(); ++it)
		{
			// count bonds and aromatic bonds
			Size num_bonds(0), num_aro(0);
			HashSet<Bond*> bonds;

			vector<Atom*>::iterator ait1 = it->begin();
			for (; ait1 != it->end(); ++ait1)
			{
				vector<Atom*>::iterator ait2(ait1);
				++ait2;
				for (; ait2 != it->end(); ++ait2)
				{
					if ((**ait1).isBoundTo(**ait2))
					{
						++num_bonds;
						Bond* const b = (**ait1).getBond(**ait2);
						bonds.insert(b);
						if (b->getOrder() == Bond::ORDER__AROMATIC)
						{
							++num_aro;
						}
					}
				}
			}

			// estimate if ring is aromatic or not
			if (float(num_aro) / float(num_bonds) >= 0.5)
			{
				for (HashSet<Bond*>::Iterator bit = bonds.begin(); bit != bonds.end(); ++bit)
				{
					(*bit)->setOrder(Bond::ORDER__AROMATIC);
				}
			}
			else
			{
				for (HashSet<Bond*>::Iterator bit = bonds.begin(); +bit; ++bit)
				{
					if ((*bit)->getOrder() == Bond::ORDER__AROMATIC)
					{
						(*bit)->setOrder(Bond::ORDER__SINGLE);
					}
				}
			}
		}
	}


	void BuildBondsProcessor::deleteOverestimatedBonds_(AtomContainer& ac)
	{
		AtomIterator ait;
		BALL_FOREACH_ATOM(ac, ait)
		{
 			if (!ait->isBound()) continue;
			
			const Size group = ait->getElement().getGroup();
			if (group != 1 && group != 17) continue;

			Bond* min_bond = 0;
			float length= Limits<float>::max();
			HashSet<Bond*> bonds;
			for (Atom::BondIterator bit = ait->beginBond(); +bit;++bit)
			{
				if (length > bit->getLength())
				{
					min_bond = &*bit;
					length = bit->getLength();
				}
				bonds.insert(&*bit);
			}

			bonds.erase(min_bond);

			for (HashSet<Bond*>::ConstIterator it=bonds.begin(); +it; ++it)
			{
				(*it)->destroy();
			}
		}
	}
	
	bool BuildBondsProcessor::getMinBondLength_(float& length, Size an1, Size an2)
	{
		if (min_bond_lengths_.has(an1) && min_bond_lengths_[an1].has(an2))
		{
			length = min_bond_lengths_[an1][an2];
			return true;
		}
		
		if (min_bond_lengths_.has(an2) && min_bond_lengths_[an2].has(an1))
		{
			length = min_bond_lengths_[an2][an1];
			return true;
		}
		
		return false;
	}

	bool BuildBondsProcessor::getMaxBondLength_(float& length, Size an1, Size an2)
	{
		if (max_bond_lengths_.has(an1) && max_bond_lengths_[an1].has(an2))
		{
			length = max_bond_lengths_[an1][an2];
			return true;
		}

		if (max_bond_lengths_.has(an2) && max_bond_lengths_[an2].has(an1))
		{
			length = max_bond_lengths_[an2][an1];
			return true;
		}

		return false;
	}

	Bond::BondOrder BuildBondsProcessor::getNearestBondOrder_(float length, Size elem1, Size elem2)
	{
		Size e1, e2;
		if (bond_lengths_.has(elem1) && bond_lengths_[elem1].has(elem2))
		{
			e1 = elem1;
			e2 = elem2;
		}
		else
		{
			if (bond_lengths_.has(elem2) && bond_lengths_[elem2].has(elem1))
			{
				e1 = elem2;
				e2 = elem1;
			}
			else
			{
				// should never occur, otherwise the input paramters file is corrupted
				Log.error() << "cannot find right bond order: " << elem1 << " " << elem2 << " " << length << endl;
				return Bond::ORDER__UNKNOWN;
			}
		}
	
		Bond::BondOrder order = Bond::ORDER__UNKNOWN;
		float min_dist(Limits<float>::max());
		HashMap<Bond::BondOrder, float> bonds = bond_lengths_[e1][e2];
		HashMap<Bond::BondOrder, float>::ConstIterator it=bonds.begin();
		for (; +it; ++it)
		{
			if (min_dist > Maths::abs(it->second-length))
			{
				min_dist = Maths::abs(it->second-length);
				order = it->first;
			}
		}
		return order;
	}

	void BuildBondsProcessor::setBondLengths(const String& file_name) throw(Exception::FileNotFound)
	{
		bond_lengths_.clear();
		min_bond_lengths_.clear();
		max_bond_lengths_.clear();
		max_length_ = 0;
		readBondLengthsFromFile_(file_name);
	}

	void BuildBondsProcessor::readBondLengthsFromFile_(const String& file_name) throw(Exception::FileNotFound)
	{

		options.setDefault(BuildBondsProcessor::Option::BONDLENGTHS_FILENAME,
													 BuildBondsProcessor::Default::BONDLENGTHS_FILENAME);
		options.setDefaultBool(BuildBondsProcessor::Option::REESTIMATE_BONDORDERS_RINGS,
													 BuildBondsProcessor::Default::REESTIMATE_BONDORDERS_RINGS);
		options.setDefaultBool(BuildBondsProcessor::Option::DELETE_EXISTING_BONDS,
													 BuildBondsProcessor::Default::DELETE_EXISTING_BONDS);
		options.setDefaultBool(BuildBondsProcessor::Option::DELETE_OVERESTIMATED_BONDS,
													 BuildBondsProcessor::Default::DELETE_OVERESTIMATED_BONDS);
		
		// test file or set default file
		String filename(file_name);
		if (file_name == "")
		{
			filename = String(options.get(BuildBondsProcessor::Option::BONDLENGTHS_FILENAME));
		}

		Path path;
		String filepath = path.find(filename);
		if (filepath == "")
		{
			throw Exception::FileNotFound(__FILE__, __LINE__, filename);
		}
	
		// read the resource file
		ResourceFile * resource_db = new ResourceFile(filepath);
		if (!resource_db->isValid())
		{
			delete resource_db;
			throw Exception::FileNotFound(__FILE__, __LINE__, filepath);
		}

		// put content of file into ResourceEntry
		ResourceEntry * tree = new ResourceEntry();
		tree->mergeChildrenOf(resource_db->getRoot());
		resource_db->close();
		delete resource_db;

		// get content and test validity of root entry
		ResourceEntry* entry1;
		entry1 = tree->getRoot().findChild("BondLengthsDB");

		if (entry1 == 0)
		{
			delete tree;
			throw Exception::ParseError(__FILE__, __LINE__, filename, "BondLengthsDB");			
		}

		// bond orders and its names in the resource file
		HashMap<String, Bond::BondOrder> name_to_order;
		name_to_order["single"] = Bond::ORDER__SINGLE;
		name_to_order["double"] = Bond::ORDER__DOUBLE;
		name_to_order["triple"] = Bond::ORDER__TRIPLE;
		name_to_order["aromatic"] = Bond::ORDER__AROMATIC;

		// iterator over the whole BondLengthDB section of the tree
		ResourceEntry::Iterator it(tree->getEntry("/BondLengthsDB")->begin());
		++it;
		for (; +it; ++it)
		{
			// test if we are in the bond partner depth
			if (it->getDepth() != 2) continue;
		
			const Size an1 = it->getKey().toUnsignedInt();
			entry1 = it->getEntry("Partners");

			// iterate over all bond partners
			ResourceEntry::Iterator entry_it;
			for (entry_it = ++entry1->begin(); +entry_it; ++entry_it)
			{
				// test if we are in the bond length section depth
				if (entry_it->getDepth() != entry1->getDepth() + 1) continue;
			
				const Size an2 = entry_it->getKey().toUnsignedInt();
				ResourceEntry* entry2 = entry_it->getEntry("BondLengths");

				// iterator over all bond lengths and orders
				ResourceEntry::Iterator bond_it;
				for (bond_it = ++entry2->begin(); +bond_it; ++bond_it)
				{
					if (bond_it->getDepth() != entry2->getDepth() + 1) continue;
				
					// now fill the bond_lengths, min and, max structs
					const String& key = bond_it->getKey();
					Bond::BondOrder order = Bond::ORDER__UNKNOWN;
					const float value = bond_it->getValue().toFloat();

					if (name_to_order.has(key))
					{
						order = name_to_order[key];
						bond_lengths_[an1][an2][order] = value;
						continue;
					}
					
					if (key == "max")
					{
						if (value > max_length_) max_length_ = value;
						max_bond_lengths_[an1][an2] = value * value;

						continue;
					}

					if (key == "min")
					{
						if (value == 0) 
						{
							min_bond_lengths_[an1][an2] = LONG_MAX;
							continue;
						}

						min_bond_lengths_[an1][an2] = value * value;
					}
					else
					{
						Log.error() << "Read unknown entry: " << key << endl;
					}
				}
			}
		}

		delete tree;
	}
	
	
} // namespace BALL
