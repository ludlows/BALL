// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// $Id: rotamerLibrary.C,v 1.29.4.3 2005/07/19 14:04:57 oliver Exp $
//

#include <BALL/STRUCTURE/rotamerLibrary.h>

#include <BALL/DATATYPE/regularExpression.h>
#include <BALL/SYSTEM/file.h>
#include <BALL/SYSTEM/path.h>
#include <BALL/STRUCTURE/geometricTransformations.h>
#include <BALL/STRUCTURE/geometricProperties.h>
#include <BALL/STRUCTURE/structureMapper.h>
#include <BALL/MATHS/matrix44.h>
#include <BALL/KERNEL/bondIterator.h>
#include <BALL/KERNEL/bond.h>

using namespace std;

namespace BALL
{
	Rotamer::Rotamer()
		: P(0.0),
			chi1(0.0),
			chi2(0.0),
			chi3(0.0),
			chi4(0.0)
	{
	}

	Rotamer::Rotamer(const Rotamer& rotamer)
		:	P(rotamer.P),
			chi1(rotamer.chi1),
			chi2(rotamer.chi2),
			chi3(rotamer.chi3),
			chi4(rotamer.chi4)
	{
	}

	Rotamer::Rotamer(float new_P, float new_chi1, float new_chi2, float new_chi3, float new_chi4)
		: P(new_P),
			chi1(new_chi1),
			chi2(new_chi2),
			chi3(new_chi3),
			chi4(new_chi4)
	{
	}

	RotamerLibrary::RotamerLibrary()
		:	variants_(),
			valid_(false)
	{
	}

	RotamerLibrary::RotamerLibrary(const String& filename, const FragmentDB& fragment_db)
		:	variants_(),
			valid_(false)
	{
		valid_ = readSQWRLLibraryFile(filename, fragment_db);
	}

	RotamerLibrary::RotamerLibrary(const RotamerLibrary& library)
		:	variants_(library.variants_),
			valid_(library.valid_)
	{
	}

	RotamerLibrary& RotamerLibrary::operator = (const RotamerLibrary& rhs)
	{
		// Avoid self assignment...
		if (&rhs != this)
		{ // ...and copy all attributes.
			variants_ = rhs.variants_;
			valid_ = rhs.valid_;
		}

		return *this;
	}

	RotamerLibrary::~RotamerLibrary()
	{
		valid_ = false;
		variants_.clear();
	}

	bool RotamerLibrary::isValid() const
	{
		return valid_;
	}
	
	ResidueRotamerSet* RotamerLibrary::getRotamerSet(const String& name)
	{
		vector<ResidueRotamerSet>::iterator it = variants_.begin();

		for(; it != variants_.end(); ++it)
		{
			if(name == it->getName())
			{
				return &(*it);  
			}
		}

		// the for-loop search was not successful 
		return  0; 
	}

	Size RotamerLibrary::getNumberOfVariants() const 
	{
		return (Size)variants_.size();
	}

	Size RotamerLibrary::getNumberOfRotamers() const 
	{
		Size number = 0;
		vector<ResidueRotamerSet>::const_iterator it = variants_.begin();
		for (; it != variants_.end(); ++it)
		{
			number += it->getNumberOfRotamers();
		}

		return number;
	}

	bool RotamerLibrary::readSQWRLLibraryFile(const String& filename, const FragmentDB& fragment_db)
	{
		// clear old stuff
		variants_.clear();
		valid_ = false;

		// try to find the rotamer library file
		Path path;
		String absolute_filename = path.find(filename);

		File f(absolute_filename);
		
		// abort if the file could not be opened
		if (!f)
		{	
			throw Exception::FileNotFound(__FILE__, __LINE__, filename);
		}

		// read the file into a vector of Strings to reparse
		// it faster. Since the SQWRL format is a pain in the ass
		// we have to identify relevant lines by a regular expression
		HashSet<String> amino_acids;
		vector<String>	lines;
		String line;
		RegularExpression regexp("[A-Z][A-Z][A-Z] [0-9] [0-9] [0-9] [0-9] *[0-9]* *[0-9]* *[0-9\\.]*");
		while (f.good())
		{
			line.getline(f);
			if (regexp.match(line))
			{	
				lines.push_back(line);
				if (!amino_acids.has(line.getField(0)))
				{
					amino_acids.insert(line.getField(0));
				}
			}
		}
		
		// build list of variant names from template DB
		// for each amino acid from the library file
		vector<String> variant_names;
		
		
		HashSet<String>::ConstIterator set_it = amino_acids.begin();
		list<String>::const_iterator list_it;
		String s[18];
		String name, variant_name;
		for (; set_it != amino_acids.end() ; set_it++)
		{
			name = *set_it;
			list<String> variant_names = fragment_db.getVariantNames(name);
			for (list_it = variant_names.begin(); list_it != variant_names.end(); ++list_it)
			{
				bool created = false;
				Size number_of_torsions = 0;
				for (Size i = 0; i < lines.size(); i++)
				{
					line = lines[i];
					if (line(0, 3) == name)
					{
						Size number_of_fields = lines[i].split(s, 18);
						float P = line.getField(7).toFloat() / 100.0;
						Angle chi1, chi2, chi3, chi4;
						chi1.set(line.getField(11).toFloat(), false);
						number_of_torsions = 1;
						if (number_of_fields > 13)
						{
							chi2.set(line.getField(13).toFloat(), false);
							number_of_torsions = 2;
						}
						if (number_of_fields > 15)
						{
							chi3.set(line.getField(15).toFloat(), false);
							number_of_torsions = 3;
						}
						if (number_of_fields > 17)
						{
							chi4.set(line.getField(17).toFloat(), false);
							number_of_torsions = 4;
						}
						
						if (!created)
						{
							ResidueRotamerSet rs(*fragment_db.getResidue(*list_it), number_of_torsions);
							rs.addRotamer(Rotamer(P, chi1, chi2, chi3, chi4));
							variants_.push_back(rs);
							created = true;
						}
						else 
						{
							variants_[variants_.size() - 1].addRotamer(Rotamer(P, chi1, chi2, chi3, chi4));
						}
					}
				}
			}
		}

		valid_ = true;
		return true;
	}

	// default constructor for ResidueRotamerSet
	ResidueRotamerSet::ResidueRotamerSet()
		: valid_(false),
			name_(""),
			side_chain_(),
			atom_name_map_(),
			rotamers_(),
			movable_atoms_chi1_(),
			movable_atoms_chi2_(),
			movable_atoms_chi3_(),
			movable_atoms_chi4_(),
			number_of_torsions_(0)
	{
	}
        
	// constructor  for ResidueRotamerSet
	ResidueRotamerSet::ResidueRotamerSet(const Residue& residue, Size number_of_torsions)
	{
		// Test if there are some atoms
		if ((residue.countAtoms()) == 0)
		{
			Log.error()  << " ResidueRotamerSet: Residue does not contain any atoms " << std::endl;
			Log.error()  <<  "ResidueRotamerSet: Maybe the residue is not contained in the fragment database..." << std::endl;  

			valid_ = false;
			return;
		}

		name_ = residue.getName();
		side_chain_ = residue;
		number_of_torsions_ = number_of_torsions;		

		// build hash map and store the original atom coordinates
		AtomIterator atom_it = side_chain_.beginAtom();
		for (; atom_it != side_chain_.endAtom(); ++atom_it)
		{
			atom_name_map_.insert(atom_it->getName(),&(*atom_it));
			original_coordinates_.push_back(atom_it->getPosition());
		}

		// Create a hash set to store the atoms that mey not be included in the torsion.
		// These are e.g. the backbone atoms or already assigned torsions (problem case: PRO!)
		HashSet<String> assigned_atoms;

		// identify anchor_atoms: those are the atoms used to "anchor" the first torsion to.
		if (atom_name_map_.has("CA") && atom_name_map_.has("C") && atom_name_map_.has("N")) 
		{
			anchor_atoms_[0] = atom_name_map_["CA"];
			anchor_atoms_[1] = atom_name_map_["C"];
			anchor_atoms_[2] = atom_name_map_["N"];
			assigned_atoms.insert("CA");
			assigned_atoms.insert("C");
			assigned_atoms.insert("N");
		} 
		else 
		{
			Log.error() << " ResidueRotamerSet: An anchor atom is missing. " << endl;
			valid_ = false;
			return;
		}

		// Identify atoms involved in the torsion for chi1 - chi4.
		// and store them in movable_atoms_chi1-4 (first four entries each).

		// Determine the atoms for chi1 and the set of movable atoms
		if (atom_name_map_.has("N") && atom_name_map_.has("CA") && atom_name_map_.has("CB") 
				&& (number_of_torsions_ > 0))
		{
			movable_atoms_chi1_.push_back("N");
			movable_atoms_chi1_.push_back("CA");
			movable_atoms_chi1_.push_back("CB");

			if (atom_name_map_.has("CG"))  movable_atoms_chi1_.push_back("CG");
			if (atom_name_map_.has("CG1")) movable_atoms_chi1_.push_back("CG1");
			if (atom_name_map_.has("OG"))  movable_atoms_chi1_.push_back("OG");
			if (atom_name_map_.has("OG1")) movable_atoms_chi1_.push_back("OG1"); 
			if (atom_name_map_.has("NG"))  movable_atoms_chi1_.push_back("NG");
			if (atom_name_map_.has("NG1")) movable_atoms_chi1_.push_back("NG1"); 
			if (atom_name_map_.has("SG"))  movable_atoms_chi1_.push_back("SG");
			if (atom_name_map_.has("SG1")) movable_atoms_chi1_.push_back("SG1"); 

			if (movable_atoms_chi1_.size() > 3)
			{
				addMovable_(movable_atoms_chi1_, *(atom_name_map_[movable_atoms_chi1_[2]]), assigned_atoms); 
				addMovable_(movable_atoms_chi1_, *(atom_name_map_[movable_atoms_chi1_[3]]), assigned_atoms);
			}
			else
			{
				// If we couldn't find the atoms, just remember this and forget about the torsion.
				number_of_torsions_ = 0;
			}
		} 
		else
		{
			number_of_torsions_ = 0;
		}

		// Determine the atoms for chi2 and the set of movable atoms
		if (movable_atoms_chi1_.size() >= 4 && (number_of_torsions_ > 1))
		{
			movable_atoms_chi2_.push_back(movable_atoms_chi1_[1]);
			movable_atoms_chi2_.push_back(movable_atoms_chi1_[2]);
			movable_atoms_chi2_.push_back(movable_atoms_chi1_[3]);

			if (atom_name_map_.has("CD"))  movable_atoms_chi2_.push_back("CD");
			if (atom_name_map_.has("CD1")) movable_atoms_chi2_.push_back("CD1");
			if (atom_name_map_.has("OD"))  movable_atoms_chi2_.push_back("OD");
			if (atom_name_map_.has("OD1")) movable_atoms_chi2_.push_back("OD1"); 
			if (atom_name_map_.has("ND"))  movable_atoms_chi2_.push_back("ND");
			if (atom_name_map_.has("ND1")) movable_atoms_chi2_.push_back("ND1"); 
			if (atom_name_map_.has("SD"))  movable_atoms_chi2_.push_back("SD");
			if (atom_name_map_.has("SD1")) movable_atoms_chi2_.push_back("SD1"); 

			if (movable_atoms_chi2_.size() > 3)
			{
				addMovable_(movable_atoms_chi2_, *(atom_name_map_[movable_atoms_chi2_[2]]), assigned_atoms); 
				addMovable_(movable_atoms_chi2_, *(atom_name_map_[movable_atoms_chi2_[3]]), assigned_atoms); 
			}
			else
			{
				// If we couldn't find the atoms, just remember this and forget about the torsion.
				number_of_torsions_ = 1;
			}
		}
		// Determine the atoms for chi3 and the set of movable atoms
		if (movable_atoms_chi2_.size() >= 4 && (number_of_torsions_ > 2))
		{
			movable_atoms_chi3_.push_back(movable_atoms_chi2_[1]);
			movable_atoms_chi3_.push_back(movable_atoms_chi2_[2]);
			movable_atoms_chi3_.push_back(movable_atoms_chi2_[3]);

			if (atom_name_map_.has("CE"))  movable_atoms_chi3_.push_back("CE");
			if (atom_name_map_.has("CE1")) movable_atoms_chi3_.push_back("CE1");
			if (atom_name_map_.has("OE"))  movable_atoms_chi3_.push_back("OE");
			if (atom_name_map_.has("OE1")) movable_atoms_chi3_.push_back("OE1"); 
			if (atom_name_map_.has("NE"))  movable_atoms_chi3_.push_back("NE");
			if (atom_name_map_.has("NE1")) movable_atoms_chi3_.push_back("NE1"); 
			if (atom_name_map_.has("SE"))  movable_atoms_chi3_.push_back("SE");
			if (atom_name_map_.has("SE1")) movable_atoms_chi3_.push_back("SE1"); 

			if (movable_atoms_chi3_.size() > 3)
			{
				addMovable_(movable_atoms_chi3_, *(atom_name_map_[movable_atoms_chi3_[2]]), assigned_atoms); 
				addMovable_(movable_atoms_chi3_, *(atom_name_map_[movable_atoms_chi3_[3]]), assigned_atoms); 
			}
			else
			{
				// If we couldn't find the atoms, just remember this and forget about the torsion.
				number_of_torsions_ = 2;
			}
		}
		// Determine the atoms for chi4 and the set of movable atoms
		if (movable_atoms_chi3_.size() >= 4 && (number_of_torsions > 3))
		{
			movable_atoms_chi4_.push_back(movable_atoms_chi3_[1]);
			movable_atoms_chi4_.push_back(movable_atoms_chi3_[2]);
			movable_atoms_chi4_.push_back(movable_atoms_chi3_[3]);

			if (atom_name_map_.has("CZ"))  movable_atoms_chi4_.push_back("CZ");
			if (atom_name_map_.has("OZ"))  movable_atoms_chi4_.push_back("OZ");
			if (atom_name_map_.has("NZ"))  movable_atoms_chi4_.push_back("NZ");
			if (atom_name_map_.has("SZ"))  movable_atoms_chi4_.push_back("SZ");

			if (movable_atoms_chi4_.size() > 3)
			{
				addMovable_(movable_atoms_chi4_, *(atom_name_map_[movable_atoms_chi4_[2]]), assigned_atoms); 
				addMovable_(movable_atoms_chi4_, *(atom_name_map_[movable_atoms_chi4_[3]]), assigned_atoms); 
			}
			else
			{
				// If we couldn't find the atoms, just remember this and forget about the torsion.
				number_of_torsions_ = 3;
			}
		}
		
		valid_ = true;
	}


	// Copy Constructor
	ResidueRotamerSet::ResidueRotamerSet(const ResidueRotamerSet& residue_rotamer_set)
		: valid_(residue_rotamer_set.valid_),
			name_(residue_rotamer_set.name_),
			side_chain_(residue_rotamer_set.side_chain_),
			atom_name_map_(),
			rotamers_(residue_rotamer_set.rotamers_),
			movable_atoms_chi1_(residue_rotamer_set.movable_atoms_chi1_),
			movable_atoms_chi2_(residue_rotamer_set.movable_atoms_chi2_),
			movable_atoms_chi3_(residue_rotamer_set.movable_atoms_chi3_),
			movable_atoms_chi4_(residue_rotamer_set.movable_atoms_chi4_),
			number_of_torsions_(residue_rotamer_set.number_of_torsions_),
			original_coordinates_(residue_rotamer_set.original_coordinates_)
	{
		for (AtomIterator it = side_chain_.beginAtom(); +it; ++it)
		{
			atom_name_map_.insert(it->getName(), &(*it));
		}

		// identify anchor_atoms
		if (atom_name_map_.has("CA")  && atom_name_map_.has("C") && atom_name_map_.has("N")) 
		{
			anchor_atoms_[0] = atom_name_map_["CA"];
			anchor_atoms_[1] = atom_name_map_["C"];
			anchor_atoms_[2] = atom_name_map_["N"];
		} 
		else 
		{
			Log.error() << " ResidueRotamerSet: An anchor atom is missing. " << endl;
			valid_ = false;
			return;
		}
	}

	ResidueRotamerSet::~ResidueRotamerSet()
	{
	}

	const ResidueRotamerSet& ResidueRotamerSet::operator = (const ResidueRotamerSet& residue_rotamer_set)
	{
		if (this != &residue_rotamer_set)
		{
			name_ = residue_rotamer_set.getName();
			side_chain_ = residue_rotamer_set.side_chain_;
			original_coordinates_ = residue_rotamer_set.original_coordinates_;
			atom_name_map_ = residue_rotamer_set.atom_name_map_;
			number_of_torsions_ = residue_rotamer_set.number_of_torsions_;

			for (AtomIterator it = side_chain_.beginAtom(); +it; ++it)
			{
				atom_name_map_.insert(it->getName(), &(*it));
			}

			// identify anchor_atoms
			if (atom_name_map_.has("CA")  && atom_name_map_.has("C") && atom_name_map_.has("N")) 
			{
				anchor_atoms_[0] = atom_name_map_["CA"];
				anchor_atoms_[1] = atom_name_map_["C"];
				anchor_atoms_[2] = atom_name_map_["N"];
			}
			else 
			{
				Log.error() << " ResidueRotamerSet: An anchor atom is missing. " << endl;
				valid_ = false;
				return (*this);
			}

			rotamers_ = residue_rotamer_set.rotamers_;
			movable_atoms_chi1_ = residue_rotamer_set.movable_atoms_chi1_;
			movable_atoms_chi2_ = residue_rotamer_set.movable_atoms_chi2_;
			movable_atoms_chi3_ = residue_rotamer_set.movable_atoms_chi3_;
			movable_atoms_chi4_ = residue_rotamer_set.movable_atoms_chi4_;
		}

		return (*this);
	}

	// return the number of rotamers
	Size ResidueRotamerSet::getNumberOfRotamers() const
	{
		return (Size)rotamers_.size();
	}

	// return the number of torsions
	Size ResidueRotamerSet::getNumberOfTorsions() const
	{
		return number_of_torsions_;
	}

	// set the number of torsions
	void ResidueRotamerSet::setNumberOfTorsions(Size number_of_torsions)
		throw(Exception::IndexOverflow)
	{
		if (number_of_torsions > 4)
		{
			throw Exception::IndexOverflow(__FILE__, __LINE__, number_of_torsions, 4);
		}
		number_of_torsions_ = number_of_torsions;
	}

	// adds a rotamer to a ResidueRotamerSet 
	void ResidueRotamerSet::addRotamer(const Rotamer& rotamer)
  {
		rotamers_.push_back(rotamer);
	}

	/** This method modifies the atom position of the template side chain
			and then returns a copy of that.
	*/
	Residue* ResidueRotamerSet::buildRotamer(const Rotamer& rotamer)
	{
		// restore original atom coordinates (see setRotamer)
		AtomIterator atom_it = side_chain_.beginAtom();
		for (Size index = 0; +atom_it && (index < original_coordinates_.size()); ++atom_it, ++index)
		{
			atom_it->setPosition(original_coordinates_[index]);
		}

		// Transform the residue template side_chain_ such that the torsion angles of the template
		// are set to the values stored in rotamer
		setTorsionAngle_(movable_atoms_chi1_, rotamer.chi1);
		setTorsionAngle_(movable_atoms_chi2_, rotamer.chi2);
		setTorsionAngle_(movable_atoms_chi3_, rotamer.chi3);
		setTorsionAngle_(movable_atoms_chi4_, rotamer.chi4);

		Residue* r = new Residue(side_chain_); 
		return r; 
	}

	// Get the name of the residue rotamer set
	const String& ResidueRotamerSet::getName() const 
	{
		return name_;
	}

	// Set the name of the residue rotamer set
	void ResidueRotamerSet::setName(const String& name)
	{
		name_ = name;
	}

	// Test if the residue rotamer set is valid
	bool ResidueRotamerSet::isValid() const 
	{
		return valid_;
	}

	// Return the residue of the residue rotamer set
	Residue& ResidueRotamerSet::getResidue()
	{
		return side_chain_;
	}
	
	// Return the residue of the residue rotamer set
	const Residue& ResidueRotamerSet::getResidue() const
	{
		return side_chain_;
	}
	
	void ResidueRotamerSet::resetResidue() 
	{
		// restore the original atom coordinates
		AtomIterator atom_it = side_chain_.beginAtom();
		for (Size index = 0; +atom_it && (index < original_coordinates_.size()); ++atom_it, ++index)
		{
			atom_it->setPosition(original_coordinates_[index]);	
		}
	}

	// Function for filling the vectors of movable atoms
	void ResidueRotamerSet::addMovable_(vector<String>& movable, Atom& a, const HashSet<String>& assigned_atoms)
	{
		Atom::BondIterator bond_it = a.beginBond();
		
		for (; bond_it != a.endBond(); ++bond_it)
		{
			Atom* b = (*bond_it).getPartner(a);
			String pattern = b->getName();

			if (find(movable.begin(), movable.end(), pattern) == movable.end()
					&& !assigned_atoms.has(pattern))
			{
				movable.push_back(pattern);
				addMovable_(movable, *b, assigned_atoms);
			}
		}
	}

	// Determine a transformation that sets the torsion of the 4 atoms to the given angle
	void ResidueRotamerSet::setTorsionAngle_(const vector<String>& movable, Angle angle)
	{

		// Test if there is a real torsion
		if (movable.size() < 4)
		{
			return ;
		}

		// Initiate Transformation Processor
		TransformationProcessor proc;

		// Compute transformation that moves the torsion atoms into normal position
		Vector3 a1(atom_name_map_[movable[0]]->getPosition());
		Vector3 a2(atom_name_map_[movable[1]]->getPosition());
		Vector3 a3(atom_name_map_[movable[2]]->getPosition());
		Vector3 a4(atom_name_map_[movable[3]]->getPosition());
		Vector3	v1(0.0,0.0,0.0);
		Vector3 v2(1.0,0.0,0.0);
		Vector3 v3(0.0,1.0,0.0);

		Matrix4x4 M = StructureMapper::matchPoints(a2, a3, a1, v1, v2, v3);

		// Apply this transformation to all atoms of the residue
		proc.setTransformation(M);
		side_chain_.apply(proc); 

		// Compute the rotation matrix (rotation around the x-axes)
		// that sets the torsion angle

		Matrix4x4 R;
		R.m11 = 1.0;	
		
		Vector3 proj_a4(M * a4);
		proj_a4.x = 0.0;
		
		float length_proj_a4 = proj_a4.getLength();

		if (length_proj_a4 != 0.0)
		{
			// Normalize the projection of a4
			proj_a4 /= length_proj_a4;
			
			// Compute cos and sin of the torsion angle 
			float cos_angle = cos(angle);
			float sin_angle = sin(angle);

			// Calculate rotation matrix
			if (proj_a4.y != 0.0)
			{
				R.m23 = (cos_angle * proj_a4.z / proj_a4.y - sin_angle) / (proj_a4.y + proj_a4.z * proj_a4.z / proj_a4.y); 
				R.m22 = (cos_angle - R.m23 * proj_a4.z) / proj_a4.y;	
			}
			else 
			{
				R.m23 = cos_angle / proj_a4.z;
				R.m22 = sin_angle / proj_a4.z;
			}

			R.m32 = -R.m23;
			R.m33 =  R.m22;

			// Rotate all atoms with index larger equal 3
			for (Size i = 3; i < movable.size(); i++)
			{
				(atom_name_map_[movable[i]])->setPosition((R * (atom_name_map_[movable[i]])->getPosition()));
			}
		}
							
	}

	// Transform the side chain such that the torsion angles are identical to the angles of the rotamer
	bool ResidueRotamerSet::setRotamer(Residue& residue, const Rotamer& rotamer) 
	{
		// restore the original side chain atom coordinates (necessary to reproduce rotamers
		// with the desired precision).
		// The application of many successive transformations on the same residue
		// is possible, but leads to slightly different coordinates for
		// the same rotamer at different times.
		AtomIterator atom_it = side_chain_.beginAtom();
		for (Size index = 0; +atom_it && (index < original_coordinates_.size()); ++atom_it, ++index)
		{
			atom_it->setPosition(original_coordinates_[index]);	
		}
		
		// Transform the residue template side_chain_ such that the torsion angles of the template
		// are set to the values stored in rotamer
		if (number_of_torsions_ > 0)
		{
			setTorsionAngle_(movable_atoms_chi1_, rotamer.chi1);
		}
		if (number_of_torsions_ > 1)
		{
			setTorsionAngle_(movable_atoms_chi2_, rotamer.chi2);
		}
		if (number_of_torsions_ > 2)
		{
			setTorsionAngle_(movable_atoms_chi3_, rotamer.chi3);
		}
		if (number_of_torsions_ > 3)
		{
			setTorsionAngle_(movable_atoms_chi4_, rotamer.chi4);
		}

		// Initiate vectors
		Vector3 a1;
		Vector3 a2;
		Vector3 a3;
		Vector3 z1;
		Vector3 z2;
		Vector3 z3;
		
		// Initiate counter
		Size counter = 0;

		// Search for the 3 backbone atoms that are needed for calculating the transformation
		for (atom_it = residue.beginAtom(); atom_it != residue.endAtom(); ++atom_it)
		{
			if ((*atom_it).getName() == "CB")
			{
				z1 = (*atom_it).getPosition();
				counter++;
			}
			if ((*atom_it).getName() == "CA") 
			{
				z2 = (*atom_it).getPosition();
				counter++;
			}
			if ((*atom_it).getName() == "N") 
			{
				z3 = (*atom_it).getPosition();
				counter++;
			}
		}

		if (counter < 3) 
		{
			return false;
		}

		// Check if the residue template side_chain_ contains the 3 backbone atoms for the matching
		if (atom_name_map_.has("CB") && atom_name_map_.has("CA") && atom_name_map_.has("N"))
		{
			a1 = atom_name_map_["CB"]->getPosition();
			a2 = atom_name_map_["CA"]->getPosition();
			a3 = atom_name_map_["N"]->getPosition();
		}
		else 
		{
			Log.error() << " template of " << side_chain_.getName() 
									<< " does not contain all backbone atoms " << endl;
			return false;
		}

		// Initiate transformation processor
		TransformationProcessor proc;

		// Calculate the transformation
		Matrix4x4 M = StructureMapper::matchPoints(a1,a2,a3,z1,z2,z3);
	
		// Apply the transformation to side_chain
		proc.setTransformation(M);
		side_chain_.apply(proc);

		// Change the coordinates of residue
		for (atom_it = residue.beginAtom(); +atom_it; ++atom_it)
		{
			if (atom_name_map_.has(atom_it->getName()))
			{
				for (Size i = 3; i < movable_atoms_chi1_.size(); ++i)
				{
					if (movable_atoms_chi1_[i] == atom_it->getName())
					{
						atom_it->setPosition(atom_name_map_[atom_it->getName()]->getPosition());
						break;
					}
				}
			}
			else 
			{
				Log.error() << " missing template atom " << atom_it->getName() << endl;
				return false;
			}
		}

		return true;
	}


	Rotamer ResidueRotamerSet::getRotamer(const Residue& residue)  const
	{
		Rotamer r;
		r.P = 1.0;

		const Atom* a1 = 0;
		const Atom* a2 = 0;
		const Atom* a3 = 0;
		const Atom* a4 = 0;

		AtomConstIterator atom_it;
		Size count;

		if (number_of_torsions_ > 0)
		{
			count = 0;
			for (atom_it = residue.beginAtom(); +atom_it; ++atom_it)
			{
				if (atom_it->getName() == movable_atoms_chi1_[0])
				{
					a1 = &*atom_it;
					count++;
				}
				if (atom_it->getName() == movable_atoms_chi1_[1])
				{
					a2 = &*atom_it;
					count++;
				}
				if (atom_it->getName() == movable_atoms_chi1_[2])
				{
					a3 = &*atom_it;
					count++;
				}
				if (atom_it->getName() == movable_atoms_chi1_[3])
				{
					a4 = &*atom_it;
					count++;
				}
			}
			
			if (count == 4)
			{
				r.chi1 = calculateTorsionAngle(*a1, *a2, *a3, *a4);
			}
		}

	
		if (number_of_torsions_ > 1)
		{
			count = 0;
			for (atom_it = residue.beginAtom(); +atom_it; ++atom_it)
			{
				if (atom_it->getName() == movable_atoms_chi2_[0])
				{
					a1 = &*atom_it;
					count++;
				}
				if (atom_it->getName() == movable_atoms_chi2_[1])
				{
					a2 = &*atom_it;
					count++;
				}
				if (atom_it->getName() == movable_atoms_chi2_[2])
				{
					a3 = &*atom_it;
					count++;
				}
				if (atom_it->getName() == movable_atoms_chi2_[3])
				{
					a4 = &*atom_it;
					count++;
				}
			}
			
			if (count == 4)
			{
				r.chi2 = calculateTorsionAngle(*a1, *a2, *a3, *a4);
			}
		}

	
		if (number_of_torsions_ > 2)
		{
			count = 0;
			for (atom_it = residue.beginAtom(); +atom_it; ++atom_it)
			{
				if (atom_it->getName() == movable_atoms_chi3_[0])
				{
					a1 = &*atom_it;
					count++;
				}
				if (atom_it->getName() == movable_atoms_chi3_[1])
				{
					a2 = &*atom_it;
					count++;
				}
				if (atom_it->getName() == movable_atoms_chi3_[2])
				{
					a3 = &*atom_it;
					count++;
				}
				if (atom_it->getName() == movable_atoms_chi3_[3])
				{
					a4 = &*atom_it;
					count++;
				}
			}
			
			if (count == 4)
			{
				r.chi3 = calculateTorsionAngle(*a1, *a2, *a3, *a4);
			}
		}

	
		if (number_of_torsions_ > 3)
		{
			count = 0;
			for (atom_it = residue.beginAtom(); +atom_it; ++atom_it)
			{
				if (atom_it->getName() == movable_atoms_chi4_[0])
				{
					a1 = &*atom_it;
					count++;
				}
				if (atom_it->getName() == movable_atoms_chi4_[1])
				{
					a2 = &*atom_it;
					count++;
				}
				if (atom_it->getName() == movable_atoms_chi4_[2])
				{
					a3 = &*atom_it;
					count++;
				}
				if (atom_it->getName() == movable_atoms_chi4_[3])
				{
					a4 = &*atom_it;
					count++;
				}
			}
			
			if (count == 4)
			{
				r.chi4 = calculateTorsionAngle(*a1, *a2, *a3, *a4);
			}
		}

		return r;
	}


	const Rotamer& ResidueRotamerSet::getRotamer(Position index) const
	{
		return rotamers_[index % getNumberOfRotamers()];
	}

#ifdef BALL_NO_INLINE_FUNCTIONS
#include <BALL/STRUCTURE/rotamerLibrary.iC>
#endif
 
} // namespace BALL
