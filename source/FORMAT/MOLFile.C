// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// $Id: MOLFile.C,v 1.25.4.1 2005/08/12 12:49:51 amoll Exp $
//


#include <BALL/FORMAT/MOLFile.h>
#include <BALL/KERNEL/atom.h>
#include <BALL/KERNEL/bond.h>
#include <BALL/KERNEL/molecule.h>
#include <BALL/KERNEL/system.h>
#include <BALL/KERNEL/PTE.h>
#include <BALL/KERNEL/forEach.h>

#define MOLFILE_VERSION_STRING_2 "V2000"
#define MOLFILE_VERSION_STRING_3 "V3000"

// enable/disable some debug output
#define DEBUG
#undef DEBUG

namespace BALL 
{

	// the format definition of the counts line
	const String MOLFile::counts_format_ = "%3d%3d%3d%3d%3d%3d%3d%3d%3d%3d%3d%6s\n";

	// the format definition for the atom block
	const String MOLFile::atom_format_ = "%10.4f%10.4f%10.4f %-3s%2d%3d%3d%3d%3d%3d%3d%3d%3d%3d%3d%3d\n";

	// the format definition for the bond block
	const String MOLFile::bond_format_ = "%3d%3d%3d%3d%3d%3d%3d\n";

	//
  const String MOLFile::Property::ATOM_MASS_DIFFERENCE = "MOLFile::MASS_DIFFERENCE";
	//
  const String MOLFile::Property::ATOM_HYDROGEN_COUNT = "MOLFile::HYDROGEN_COUNT";
	//
  const String MOLFile::Property::ATOM_STEREO_CARE_BOX = "MOLFile::STEREO_CARE_BOX";
  //
  const String MOLFile::Property::ATOM_VALENCE = "MOLFile::ATOM_VALENCE";
  //
  const String MOLFile::Property::ATOM_H0_DESIGNATOR = "MOLFile::H0_DESIGNATOR";
  //
	const String MOLFile::Property::ATOM_REACTION_COMPONENT_TYPE = "MOLFile::REACTION_COMPONENT_TYPE";
	//
	const String MOLFile::Property::ATOM_REACTION_COMPONENT_NUMBER = "MOLFile::REACTION_COMPONENT_NUMBER";
	//
	const String MOLFile::Property::ATOM_INVERSION_RETENTION = "MOLFile::INVERSION_RETENTION";
	//
	const String MOLFile::Property::ATOM_EXACT_CHANGE = "MOLFile::ATOM_EXACT_CHANGE";

	//
	const String MOLFile::Property::BOND_STEREO = "MOLFile::STEREO";
	//
	const String MOLFile::Property::BOND_TOPOLOGY = "MOLFile::TOPOLOGY";
	//
	const String MOLFile::Property::BOND_REACTING_CENTER_STATUS = "MOLFile::REACTING_CENTER_STATUS";


	MOLFile::MOLFile()
		throw()
		:	GenericMolFile()
	{
	}

	MOLFile::MOLFile(const String& name, File::OpenMode open_mode)
		throw(Exception::FileNotFound)
		: GenericMolFile()
	{
		GenericMolFile::open(name, open_mode);
	}

	MOLFile::MOLFile(const MOLFile& file)
		throw(Exception::FileNotFound)
		: GenericMolFile(file)
	{
	}

	MOLFile::~MOLFile()
		throw()
	{
	}

	bool MOLFile::write(const Molecule& molecule)
		throw(File::CannotWrite)
	{
		if (!isOpen() || getOpenMode() != std::ios::out)
		{
			throw File::CannotWrite(__FILE__, __LINE__, name_);
		}

		// write header block
		String name = molecule.getName();
		if ((name.size() > 80) || (name.has('\n')))
		{
			Log.warn() << "MOLFile::write: truncating illegal molecule name ('" << name << "')." << std::endl;
			if (name.size() > 80)
			{
				name = name(0, 80);
			}
			if (name.has('\n'))
			{
				name = name.before("\n");
			}
			if (name.hasSubstring("$$$$") || name.hasSubstring("$MDL") 
					|| name.hasSubstring("$RXN") || name.hasSubstring("$RDFILE")) 
			{
				name = "--";
			}
		}

		getFileStream() << name << std::endl;
		getFileStream() << "      " << "BALL " << VersionInfo::getVersion() << std::endl;
		getFileStream() << std::endl;
		

		// write counts line
		CountsStruct counts;
		counts.number_of_atoms = molecule.countAtoms();
		counts.number_of_bonds = molecule.countBonds();
		counts.number_of_atom_lists = 0;
		counts.chiral = false;
		counts.number_of_stext_entries = 0;
		counts.number_of_reaction_components = 0;
		counts.number_of_reactants = 0;
		counts.number_of_products = 0;
		counts.number_of_intermediates = 0;
		counts.version = MOLFILE_VERSION_STRING_2;

		writeCountsLine_(counts);
		
		// write atom block
		HashMap<const Atom*, Position> atom_map;
		AtomStruct atom;
		AtomConstIterator it = molecule.beginAtom();
		Position atom_number = 1;
		for (; +it; ++it)
		{
			atom.symbol = it->getElement().getSymbol();
			atom.position = it->getPosition();
			atom.mass_difference = 0;
			switch ((int)it->getFormalCharge())
			{
				case  3: atom.charge = 1; break;
				case  2: atom.charge = 2; break;
				case  1: atom.charge = 3; break;
				case -1: atom.charge = 5; break;
				case -2: atom.charge = 6; break;
				case -3: atom.charge = 7; break;

				default:
					atom.charge = 0;
			}
			atom.parity = 0;
			atom.hydrogen_count = 0;
			atom.stereo_care_box = false;
			// Detect atypical valences ?????
			atom.valence = 0;
			atom.H0_designator = false;
			atom.reaction_component_type = 0;
			atom.reaction_component_number = 0;
			atom.number = 0;
			atom.inversion_retention = 0;
			atom.exact_change = false;

			// store the atom index in a hash map
			atom_map.insert(std::pair<const Atom*, Position>(&*it, atom_number++));

			writeAtomLine_(atom);
		}

		// write bond block
		BondStruct bond;
		Atom::BondConstIterator bond_it;
		BALL_FOREACH_BOND(molecule, it, bond_it)
		{
			// figure out which atom indices we need
			if (atom_map.has(bond_it->getFirstAtom()))
			{
				bond.first_atom = atom_map[bond_it->getFirstAtom()];
			}
			else
			{
				Log.warn() << "MOLFile::write: ignoring bond between " << bond_it->getFirstAtom()->getFullName() 
									 << " and " << bond_it->getSecondAtom()->getFullName() << std::endl;
				continue;
			}
			if (atom_map.has(bond_it->getSecondAtom()))
			{
				bond.second_atom = atom_map[bond_it->getSecondAtom()];
			}
			else
			{
				Log.warn() << "MOLFile::write: ignoring bond between " << bond_it->getFirstAtom()->getFullName() 
									 << " and " << bond_it->getSecondAtom()->getFullName() << std::endl;
				continue;
			}
				
			// translate BALL bond type to the MOL numerical value
			switch (bond_it->getOrder())
			{
				case Bond::ORDER__SINGLE:		bond.type = 1; break;
				case Bond::ORDER__DOUBLE:		bond.type = 2; break;
				case Bond::ORDER__TRIPLE:		bond.type = 3; break;
				case Bond::ORDER__AROMATIC: bond.type = 4; break;
				default:
					bond.type = 8; // type = ANY for all other bond types
			}
			
			bond.stereo = 0;
			bond.topology = 0;
			bond.reacting_center_status = 0;
			
			writeBondLine_(bond);
		}

		// write propery section
		getFileStream() << "M  END" << std::endl;
		
		return true;
	}

	bool MOLFile::write(const System& system)
		throw(File::CannotWrite)
	{
		MoleculeConstIterator mol = system.beginMolecule();
		if (!write(*mol)) return false;
		mol++;
		if (mol != system.endMolecule())
		{
			Log.warn() << "MOLFile::write: found more than one molecule in system while writing -- all molecules after the first one are ignored!" << std::endl;
		}
		return true;
	}

	Molecule* MOLFile::readCTAB_(vector<Atom*>& atom_map)
		throw(Exception::ParseError)
	{

		#ifdef DEBUG
			Log.info() << "entering MOLFile::readCTAB_(current line = " << getLineNumber() << ")" << std::endl;
		#endif
		// read the counts line
		CountsStruct counts;
		if (!readCountsLine_(counts))
		{
			throw Exception::ParseError(__FILE__, __LINE__, 
					String("'") + getLine() + "' (line " + String(getLineNumber()) + " of '" + getName() + "')",
																				String("Cannot parse counts line: "));

		}
		#ifdef DEBUG
			Log.info() << "Counts line: " << counts.number_of_atoms 
				<< " atoms and " << counts.number_of_bonds << " bonds." << std::endl;
		#endif

		// resize the array to the number of atoms
		atom_map.resize(counts.number_of_atoms);
		Molecule* molecule = new Molecule;
		try
		{
			// read the atom block
			AtomStruct atom_struct;
			for (Position i = 0; i < counts.number_of_atoms; i++)
			{
				if (!readAtomLine_(atom_struct))
				{
					throw Exception::ParseError(__FILE__, __LINE__, 
							String("'") + getLine() + "' (line " + String(getLineNumber()) + " of '" + getName() + "')",
																						String("Cannot parse atom line: "));

				}

				// create the atom
				Atom* atom = new Atom(PTE[atom_struct.symbol.trim()], 
															atom_struct.symbol.trim(), // name
															"", // emtpy type name
															Atom::UNKNOWN_TYPE,
															atom_struct.position,
															Vector3(0.0, 0.0, 0.0), // velocity undefined
															Vector3(0.0, 0.0, 0.0), // force undefined
															0.0, // charge reset below
															0.0); // radius undefined
				// insert the atom into the molecule and store its pointer in 
				// a map array
				molecule->append(*atom);
				atom_map[i] = atom;

				// assign the charge (overridden by M CHG lines below)
				switch (atom_struct.charge)
				{
					case 1:	atom->setCharge(3.0); atom->setFormalCharge(3); break;
					case 2:	atom->setCharge(2.0); atom->setFormalCharge(2); break;
					case 3:	atom->setCharge(1.0); atom->setFormalCharge(1); break;
					// doublet -- how to handle this?
					case 4:	atom->setCharge(0.0); atom->setFormalCharge(0); break;
					case 5: atom->setCharge(-1.0); atom->setFormalCharge(-1); break;
					case 6: atom->setCharge(-2.0); atom->setFormalCharge(-2); break;
					case 7: atom->setCharge(-3.0); atom->setFormalCharge(-3); break;

					case 0: // charge 0 or other
					default:
						atom->setCharge(0.0);
						atom->setFormalCharge(0);
						break;
				}
				
				// store the remaining information as named properties
				atom->setProperty(Property::ATOM_MASS_DIFFERENCE, atom_struct.mass_difference);
				atom->setProperty(Property::ATOM_HYDROGEN_COUNT, atom_struct.hydrogen_count);
				atom->setProperty(Property::ATOM_STEREO_CARE_BOX, atom_struct.stereo_care_box);
				atom->setProperty(Property::ATOM_VALENCE, atom_struct.valence);
				atom->setProperty(Property::ATOM_H0_DESIGNATOR, atom_struct.H0_designator);
				atom->setProperty(Property::ATOM_REACTION_COMPONENT_TYPE, atom_struct.reaction_component_type);
				atom->setProperty(Property::ATOM_REACTION_COMPONENT_NUMBER, atom_struct.reaction_component_number);
				atom->setProperty(Property::ATOM_INVERSION_RETENTION, atom_struct.inversion_retention);
				atom->setProperty(Property::ATOM_EXACT_CHANGE, atom_struct.exact_change);
			}

			// read the bond block
			BondStruct bond_struct;
			for (Position i = 0; i < counts.number_of_bonds; i++)
			{
				if (!readBondLine_(bond_struct)) return false;

				// ensure the atoms referenced do exist
				if ((bond_struct.first_atom < 1) || (bond_struct.first_atom > counts.number_of_atoms))
				{
					throw Exception::ParseError(__FILE__, __LINE__, String("'") + getLine() + "' (line " + String(getLineNumber()) + " of '" + getName() + "')",
																			String("Referencing undefined atom number: ") + String(bond_struct.first_atom));
				}
				if ((bond_struct.second_atom < 1) || (bond_struct.second_atom > counts.number_of_atoms))
				{
					throw Exception::ParseError(__FILE__, __LINE__, String("'") + getLine() + "' (line " + String(getLineNumber()) + " of '" + getName() + "')",
																			String("Referencing undefined atom number: ") + String(bond_struct.second_atom));
				}
				
				// If the bond should join the atom with itself, we safely ignore this.
				// There are some SD files out there showing this specific problem...
				if (bond_struct.first_atom == bond_struct.second_atom)
				{
					continue;
				}

				// create the bond
				Bond* bond = atom_map[bond_struct.first_atom - 1]->createBond(*(atom_map[bond_struct.second_atom - 1]));
				if (bond == 0)
				{
					throw Exception::ParseError(__FILE__, __LINE__, 
							String("'") + getLine() + "' (line " + String(getLineNumber()) + " of '" + getName() + "')",
							String("Cannot create bond: ") + String(bond_struct.second_atom) + "-" + String(bond_struct.first_atom)) ;
				}
				
				// translate the bond type
				switch (bond_struct.type)
				{
					case 1: // single
						bond->setOrder(Bond::ORDER__SINGLE); break;
					
					case 2: // double	
						bond->setOrder(Bond::ORDER__DOUBLE); break;

					case 3: // 
						bond->setOrder(Bond::ORDER__TRIPLE); break;

					case 4: // 
						bond->setOrder(Bond::ORDER__AROMATIC); break;

					case 5: // single or double
					case 6: // single or aromatic
					case 7: // double or aromatic
					case 8: // any
						break; 
					
					default:
						throw Exception::ParseError(__FILE__, __LINE__, String("'") + getLine() + "' (line " + String(getLineNumber()) + " of '" + getName() + "')",
																				String("Illegal bond type: ") + String(bond_struct.type));
				}

				// store remaining stuff as named properties	
				bond->setProperty(Property::BOND_STEREO, bond_struct.stereo);
				bond->setProperty(Property::BOND_TOPOLOGY, bond_struct.topology);
				bond->setProperty(Property::BOND_REACTING_CENTER_STATUS, bond_struct.reacting_center_status);
			}

			// skip the atom list block
			for (Position i = 0; i < counts.number_of_atom_lists; i++)
			{
				readLine();
			}

			// skip the stext block
			for (Position i = 0; i < counts.number_of_stext_entries; i++)
			{
				readLine();
				readLine();
			}

			// read the properties block
			while (!startsWith("M  END") && readLine() && good() && startsWith("M "))
			{
				// delete the "M  " part of the line -- the next three letters are the tag.
				String tag(getLine().getSubstring(3, 3));
				tag.trim();

				if (tag == "END")
				{
					// end of section, abort the wile loop
					break;
				}
				if (tag == "A") // atom alias
				{
					//????
				}
				else if (tag == "V") // atom value
				{
				}
				else if (tag == "G") // group
				{
				}
				else if (tag == "CHG") // charge
				{
				}
				else if (tag == "RAD") // radical
				{
				}
				else if (tag == "ISO") // radical
				{
				}
				else if (tag == "RBC") // ring bond count
				{
				}
				else if (tag == "SUB") // substitution count
				{
				}
				else if (tag == "UNS") // unsaturated atom
				{
				}
				else if (tag == "LIN") // link atom
				{
				}
				else if (tag == "ALS") // atom list
				{
				}
				else if (tag == "APO") // attachment point
				{
				}
				else if (tag == "AAL") // atom attachment order
				{
				}
				else if (tag == "RGP") // Rgroup label location
				{
				}
				else if (tag == "LOG") // Rgroup logic, unstisfied sites, Range of occurrence
				{
				}
				else if (tag == "STY") // Sgroup type
				{
				}
				else if (tag == "SST") // Sgroup subtype
				{
				}
				else if (tag == "SLB") // Sgroup labels
				{
				}
				else if (tag == "SCN") // Sgroup connectivity
				{
				}
				else if (tag == "") // Sgroup connectivity
				{
				}
				else
				{
					Log.warn() << "MOLFile::readCTAB_: ignoring property entry " << tag << std::endl;
				}
			}
		}
		catch (Exception::GeneralException& e)
		{
			#ifdef DEBUG
				Log.info() << "MOLFile::readCTAB_: caught exception while parsing line " << getLineNumber() << ": " << e << std::endl;
			#endif
			// clean up: delete all atoms we just constructed
			delete molecule;
			molecule = 0;
			
			throw e;
		}
		
		#ifdef DEBUG
			Log.info() << "MOLFile::readCTAB_ = " << (void*)molecule << std::endl;
		#endif

		return molecule;
	}

	bool MOLFile::read(System& system)
		throw(Exception::ParseError)
	{
		// read the molecule
		Molecule* molecule = 0;
		try
		{
			molecule = read();
		}
		catch (Exception::ParseError& e)
		{
			#ifdef DEBUG
				Log.info() << "MOLFile::read(System&): caught exception while parsing line " << getLineNumber() << ": " << e << std::endl;
			#endif

			// clean up and rethrow
			delete molecule;
			throw e;
		}

		if (molecule == 0) 
		{
			return false;
		}

		// add the molecule to the system
		system.append(*molecule);

		return true;
	}

	Molecule* MOLFile::read()
		throw(Exception::ParseError)
	{
		// read the header block: first line == name, third line = comment, second line ignored
		bool ok = readLine();
		if (!ok || !good())
		{
			// end of file encountered or not open -- nothing there to read...
			return 0;
		}

		
		String name = getLine();
		ok &= readLine();
		ok &= readLine();
		String comment = getLine();
		if (!ok)
		{
			throw Exception::ParseError(__FILE__, __LINE__, String("'") + getLine() + "' (line " + String(getLineNumber()) + " of '" + getName() + "')",
																	"Unable to read header block");
		}
		static vector<Atom*> atom_map;
		Molecule* mol = readCTAB_(atom_map);
		if (mol) mol->setName(name);

		return mol;
	}
	
	bool MOLFile::readCountsLine_(MOLFile::CountsStruct& counts)
	{
		// read the next line
		readLine();

		bool ok = true;

		// parse the line according to the Counts format
		counts.number_of_atoms = 0;
		ok &= parseColumnFormat("%3d", 0, 3, (void*)&counts.number_of_atoms);

		counts.number_of_bonds = 0;
		ok &= parseColumnFormat("%3d", 3, 3, (void*)&counts.number_of_bonds);

		// We do not assume that anything except the number of atoms
		// and bonds (first two fields) *has* to be specified.
		Size len = getLine().size();
		counts.number_of_atom_lists = 0;
		counts.chiral = false;
		counts.number_of_stext_entries = 0;
		counts.number_of_reaction_components = 0;
		counts.number_of_reactants = 0;
		counts.number_of_products = 0;
		counts.number_of_intermediates = 0;

		if (len >= 9) parseColumnFormat("%3d", 6, 3, (void*)&counts.number_of_atom_lists);
		if (len >= 15)
		{
			Size chiral = 0;
			parseColumnFormat("%3d", 12, 3, (void*)&chiral);
			counts.chiral = (chiral == 0);
		}		
		if (len >= 18) parseColumnFormat("%3d", 15, 3, (void*)&counts.number_of_stext_entries);
		if (len >= 21) parseColumnFormat("%3d", 18, 3, (void*)&counts.number_of_reaction_components);
		if (len >= 24) parseColumnFormat("%3d", 21, 3, (void*)&counts.number_of_reactants);
		if (len >= 27) parseColumnFormat("%3d", 24, 3, (void*)&counts.number_of_products);
		if (len >= 30) parseColumnFormat("%3d", 27, 3, (void*)&counts.number_of_intermediates);

		return ok;
	}

	void MOLFile::writeCountsLine_(const MOLFile::CountsStruct& counts)
	{
		static char buf[81];
		sprintf(buf, counts_format_.c_str(), 
						counts.number_of_atoms,
						counts.number_of_bonds,
						counts.number_of_atom_lists,
						999, // field is obsolete
						(int)counts.chiral,
						counts.number_of_stext_entries,
						counts.number_of_reaction_components,
						counts.number_of_reactants,
						counts.number_of_products,
						counts.number_of_intermediates,
						999, // number of additional properties -- no longer supported (default = 999)
						counts.version.c_str());
		getFileStream() << buf;
	}

	bool MOLFile::readAtomLine_(MOLFile::AtomStruct& atom)
	{
		// read the next line
		readLine();

		bool ok = true;

		// parse the line according to the atom block format:
		// "%10.4f%10.4f%10.4f %3s%2d%3d%3d%3d%3d%3d%3d%3d%3d%3d%3d%3d\n"
		atom.position.x = 0.0;
		ok &= parseColumnFormat("%f", 0, 10, (void*)&(atom.position.x));
		
		atom.position.y = 0.0;
		ok &= parseColumnFormat("%f", 10, 10, (void*)&(atom.position.y));

		atom.position.z = 0.0;
		ok &= parseColumnFormat("%f", 20, 10, (void*)&(atom.position.z));
		
		atom.symbol = "?";
		atom.symbol = String(getLine(), 31, 3);

		atom.mass_difference = 0;
		ok &= parseColumnFormat("%d", 34, 2, (void*)&atom.mass_difference);

		atom.charge = 0;
		ok &= parseColumnFormat("%d", 36, 3, (void*)&atom.charge);

		atom.parity = 0;
		ok &= parseColumnFormat("%d", 39, 3, (void*)&atom.parity);

		atom.hydrogen_count = 0;
		ok &= parseColumnFormat("%d", 42, 3, (void*)&atom.hydrogen_count);

		atom.stereo_care_box = 0;
		ok &= parseColumnFormat("%d", 45, 3, (void*)&atom.stereo_care_box);

		atom.valence = 0;
		atom.reaction_component_type = 0;
		atom.reaction_component_number = 0;
		atom.H0_designator = 0;
		atom.number = 0;
		atom.inversion_retention = 0;
		atom.exact_change = 0;
		Size len = getLine().size();
		if (len >= 51) parseColumnFormat("%d", 48, 3, (void*)&atom.valence);
		if (len >= 54) parseColumnFormat("%d", 51, 3, (void*)&atom.H0_designator);
		if (len >= 57) parseColumnFormat("%d", 54, 3, (void*)&atom.reaction_component_type);
		if (len >= 60) parseColumnFormat("%d", 57, 3, (void*)&atom.reaction_component_number);
		if (len >= 63) ok &= parseColumnFormat("%d", 60, 3, (void*)&atom.number);
		if (len >= 66) ok &= parseColumnFormat("%d", 63, 3, (void*)&atom.inversion_retention);
		if (len >= 69) ok &= parseColumnFormat("%d", 66, 3, (void*)&atom.exact_change);

		return ok;
	}

	void MOLFile::writeAtomLine_(const MOLFile::AtomStruct& atom)
	{
		static char buf[80];
		sprintf(buf,  atom_format_.c_str(), 
						atom.position.x,
						atom.position.y,
						atom.position.z,
						atom.symbol.c_str(),
						atom.mass_difference,
						atom.charge,
						atom.parity,
						atom.hydrogen_count,
						(atom.stereo_care_box ? 1 : 0),
						atom.valence,
						(atom.H0_designator ? 1 : 0),
						atom.reaction_component_type,
						atom.reaction_component_number,
						atom.number,
						atom.inversion_retention,
						atom.exact_change);
		getFileStream() << buf;		
	}

	bool MOLFile::readBondLine_(MOLFile::BondStruct& bond)
	{
		// read the next line
		readLine();

		bool ok = true;

		// parse the line according to the bond block format
		bond.first_atom = 0;
		ok &= parseColumnFormat("%3d", 0, 3, (void*)&bond.first_atom);

		bond.second_atom = 0;
		ok &= parseColumnFormat("%3d", 3, 3, (void*)&bond.second_atom);

		bond.type = 0;
		ok &= parseColumnFormat("%3d", 6, 3, (void*)&bond.type);

		bond.stereo = 0;
		ok &= parseColumnFormat("%3d", 9, 3, (void*)&bond.stereo);

		bond.topology = 0;
		ok &= parseColumnFormat("%3d", 15, 3, (void*)&bond.topology);

		bond.reacting_center_status = 0;
		ok &= parseColumnFormat("%3d", 15, 3, (void*)&bond.reacting_center_status);

		return ok;
	}

	void MOLFile::writeBondLine_(const MOLFile::BondStruct& bond)
	{
		static char buf[80];
		sprintf(buf, bond_format_.c_str(),
						bond.first_atom,
						bond.second_atom,
						bond.type,
						bond.stereo,
						0,
						bond.topology,
						bond.reacting_center_status);
		getFileStream() << buf;
	}

	const MOLFile& MOLFile::operator = (const MOLFile& file)
		throw()
	{
		GenericMolFile::operator = (file);
		return *this;
	}

} // namespace BALL
