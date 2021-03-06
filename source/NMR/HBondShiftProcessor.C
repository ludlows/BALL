// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// $Id: HBondShiftProcessor.C,v 1.6 2003/08/26 09:18:09 oliver Exp $
//

#include <BALL/NMR/HBondShiftProcessor.h>
#include <stdio.h>

using namespace std;

namespace BALL
{
  const char* HBondShiftProcessor::PROPERTY__HBOND_SHIFT = "HBondShift";

  HBondShiftProcessor::HBondShiftProcessor()
    throw()
    : ShiftModule()
  {
  }

  HBondShiftProcessor::HBondShiftProcessor(const HBondShiftProcessor& processor)
    throw()
    : ShiftModule(processor),
      donor_list_(processor.donor_list_),
      acceptor_list_(processor.acceptor_list_),
      a_(processor.a_),
      b_(processor.b_),
      minimum_bond_length_(processor.minimum_bond_length_),
      maximum_bond_length_(processor.maximum_bond_length_)
  {
  }

  HBondShiftProcessor::~HBondShiftProcessor()
    throw()
  {
  }

  void HBondShiftProcessor::init()
    throw()
  {
		// by default, we assume the worst...
    valid_ = false;

    // if no parameters are assigned, abort immediately
    if (parameters_ == 0)
    {
      return;
    }

    // check that the parameter file contains the correct section...
    ParameterSection parameter_section;
    parameter_section.extractSection(*parameters_, "HBondEffect");

    // ...and that this section contains the correct parameters, and if so, assign them to variables
    if (!parameter_section.hasVariable("a") || !parameter_section.hasVariable("b")
	|| !parameter_section.hasVariable("minimum_bond_length") || !parameter_section.hasVariable("maximum_bond_length"))
    {
      return;
    }

    String numdum;
    numdum = parameter_section.getValue("a", "value");
    a_ = numdum.toFloat();
    numdum = parameter_section.getValue("b", "value");
    b_ = numdum.toFloat();
    numdum = parameter_section.getValue("minimum_bond_length", "value");
    minimum_bond_length_ = numdum.toFloat();
    numdum = parameter_section.getValue("maximum_bond_length", "value");
    maximum_bond_length_ = numdum.toFloat();
 
    // mark the module as initialized
    valid_ = true;
  }

  bool HBondShiftProcessor::start()
    throw()
  {
    // if the module is invalid, abort
    if (!isValid())
    {
      return false;
    }

    // clear the donor list and the acceptor list
    donor_list_.clear();
    acceptor_list_.clear();

    return true;
  }

  bool HBondShiftProcessor::finish()
    throw()
  {
    // if the module is invalid, abort
    if (!isValid())
    {
      return false;
    }

    // if there were no donors or acceptors, return immediately
    if ((donor_list_.size() == 0) || acceptor_list_.size() == 0)
    {
      return true;
    }

    // iterate over all donors
    list<Atom*>::iterator donor_it = donor_list_.begin();
    for (; donor_it != donor_list_.end(); ++donor_it)
    {
      Vector3 donor_pos = (*donor_it)->getPosition();
      float delta_HBOND = 0.0;
      Size hbond_number = 0;
      float nearest_acc = 0.0;
		 
      // iterate over all acceptors
      list<Atom*>::iterator acceptor_it = acceptor_list_.begin();
      for (; acceptor_it != acceptor_list_.end(); ++acceptor_it)
      { 

	// Test: use the distances of the hydrogens instead of the ones of the residues.
	BALL::Atom::BondConstIterator bi = (*donor_it)->beginBond();
	
	for (; +bi; ++bi)
	{
	  Atom* at = (*bi).getPartner(**donor_it);
	  if (at->getElement() == PTE[Element::H])
	    {
	      // print the partners of the assumed hydrogen bond
	      
	      // compute the euclidian distance between acceptor and donor
	      Vector3 acceptor_pos = (*acceptor_it)->getPosition();
	      Vector3 donor_h_pos = (at)->getPosition();
	      
	      Vector3 don_h = donor_pos - donor_h_pos;
	      Vector3 don_h_acc = acceptor_pos - donor_h_pos;
	      
	      Angle psi = don_h.getAngle(don_h_acc);
	      	      
	      BALL::Atom::BondConstIterator accit = (*acceptor_it)->beginBond();
	      Atom* c;
	      Vector3 acceptor_plane;
	      Angle theta;
	      Angle phi;
	      
	      while (+accit) 
	      {
		c = (*accit).getPartner(**acceptor_it);
		if (c->getElement() == PTE[Element::C])
		{
		  Vector3 c_pos = c->getPosition();
		  BALL::Atom::BondConstIterator c_partner_it = c->beginBond();
		  while (+c_partner_it)
		    {
		      Atom* c_partner = (*c_partner_it).getPartner(*c);
		      if (c_partner->getElement() != PTE[Element::H])
			{
			  Vector3 c_partner_pos = c_partner->getPosition();
			  if ((c_partner_pos != c_pos) && (donor_h_pos != acceptor_pos))
			    {
			      acceptor_plane = (acceptor_pos - c_pos) % (c_partner_pos - c_pos);
			      Vector3 hbond = (donor_h_pos - acceptor_pos);
			      if ((acceptor_plane != Vector3(0, 0, 0)) && (((c_pos - acceptor_pos) % acceptor_plane) != Vector3(0, 0, 0)))
				{
				  theta = acceptor_plane.getAngle(hbond);
				  Vector3 xaxis = (c_pos - acceptor_pos) % acceptor_plane;
				  phi = xaxis.getAngle(hbond);
				  break;
				}
			    }
			}
		      ++c_partner_it;
		    }
		  break;
		}
		++accit;
	      }
	      
	      
	      
	      // this is used to compute the shift	
	      float b_length = donor_h_pos.getDistance(acceptor_pos);
	      
	      // this is used to decide whether there is an hbond or not
	      float distance = donor_pos.getDistance(acceptor_pos);
	      
	      // this is a test! empirical guess
	      // float b_length = (2. / 3) * distance;
	      
	      // do they form a bond?
	      // like in the paper by wuethrich, no bonds are accepted: - outside a certain range
	      //							  - on the same residue
	      //							  - on neighbouring residues in the chain
	      if ((minimum_bond_length_ <= distance) && (maximum_bond_length_ >= distance) 
		  && (*acceptor_it)->getResidue() != (*donor_it)->getResidue()
		  && (abs((*acceptor_it)->getResidue()->getID().toInt() - (*donor_it)->getResidue()->getID().toInt()) > 1))
		{
		  // increase the number of found HBonds for this donor
		  cout << (*donor_it)->getResidue()->getFullName() << ":" << (*donor_it)->getResidue()->getID() << " -> ";
		  cout << (*acceptor_it)->getResidue()->getFullName() << ":" << (*acceptor_it)->getResidue()->getID() << " " << endl;
		  hbond_number++;
		  if (hbond_number > 1)
		    {
		      // only accept the nearest partner
		      if (b_length <= nearest_acc)
			{
			  // compute their contribution to the chemical shift
			  delta_HBOND = a_ * 1./pow(b_length, 3) + b_;
			  cout << "#: " << hbond_number << " length: " << b_length << " " << distance << " " << delta_HBOND << endl;
			  nearest_acc = b_length;
			  hbond_number = 1;
			  cerr << (*donor_it)->getResidue()->getFullName() << (*donor_it)->getResidue()->getID() << " " << (*donor_it)->getName() << " -> ";
			  cerr << (*acceptor_it)->getResidue()->getFullName() << (*acceptor_it)->getResidue()->getID() << " " << (*acceptor_it)->getName() << " " << b_length << " " << distance << " phi, theta, psi: " << phi << " " << theta << " " << psi << endl;
			}
		    }
		  else
		    {
		      cerr << (*donor_it)->getResidue()->getFullName() << (*donor_it)->getResidue()->getID() << " " << (*donor_it)->getName() << " -> ";
		      cerr << (*acceptor_it)->getResidue()->getFullName() << (*acceptor_it)->getResidue()->getID() << " " << (*acceptor_it)->getName() << " " << b_length << " " << distance << " phi, theta, psi: " << phi << " " << theta << " " << psi << endl;
		      delta_HBOND = a_ * 1./pow(b_length, 3) + b_;
		      cout << "#: " << hbond_number << " length: " << b_length << " " << distance << " " << delta_HBOND << endl;
		      nearest_acc = b_length;
		    }
		}
	    }
	}
	bi = (*donor_it)->beginBond();
	
      }
      BALL::Atom::BondConstIterator bi = (*donor_it)->beginBond();
      if (hbond_number != 0)
	{ 
	  delta_HBOND /= (float) hbond_number; 
	  for (; +bi; ++bi)
	    {
	      Atom* at = (*bi).getPartner(**donor_it);
	      if (at->getElement() == PTE[Element::H])
		{
		  float shift = (at)->getProperty(ShiftModule::PROPERTY__SHIFT).getFloat();
		  shift += delta_HBOND;	
		  (at)->setProperty(ShiftModule::PROPERTY__SHIFT, shift);
		  (at)->setProperty(PROPERTY__HBOND_SHIFT, delta_HBOND);
		}
	    }
	}
    } 
    return true;
  }
  
  Processor::Result HBondShiftProcessor::operator () (Composite& object)
    throw()
  {
    // here, we collect all possible acceptors and donors
    if (RTTI::isKindOf<Atom>(object))
      {
	Atom* atom_ptr = RTTI::castTo<Atom>(object);

	// ??? Do I have to use the fragment or the residue of the atom???
	const String& residue_name = atom_ptr->getFragment()->getName();
	const String& atom_name = atom_ptr->getName();

	if (   ((residue_name == "ALA") && (atom_name == "N"))
	       || ((residue_name == "ARG") && ((atom_name == "N") || (atom_name == "NH1") || (atom_name == "NH2")))
	       || ((residue_name == "ASN") && ((atom_name == "N") || (atom_name == "ND2")))
	       || ((residue_name == "ASP") && (atom_name == "N"))
	       || ((residue_name == "CYS") && (atom_name == "N"))
	       || ((residue_name == "GLN") && ((atom_name == "N") || (atom_name == "NE2")))
	       || ((residue_name == "GLU") && (atom_name == "N"))
	       || ((residue_name == "GLY") && (atom_name == "N"))
	       || ((residue_name == "HIS") && ((atom_name == "N") || (atom_name == "NE2")))
	       || ((residue_name == "ILE") && (atom_name == "N"))
	       || ((residue_name == "LEU") && (atom_name == "N"))
	       || ((residue_name == "LYS") && ((atom_name == "N") || (atom_name == "NZ")))
	       || ((residue_name == "MET") && (atom_name == "N"))
	       || ((residue_name == "PHE") && (atom_name == "N"))
	       || ((residue_name == "PRO") && (atom_name == "N"))
	       || ((residue_name == "SER") && ((atom_name == "N") || (atom_name == "OG")))
	       || ((residue_name == "THR") && ((atom_name == "N") || (atom_name == "OG1")))
	       || ((residue_name == "TRP") && (atom_name == "N"))
	       || ((residue_name == "TYR") && ((atom_name == "N") || (atom_name == "OH")))
	       || ((residue_name == "VAL") && (atom_name == "N")))
	  {
	    donor_list_.push_back(atom_ptr);
	  }
     
	if (   ((residue_name == "ALA") && (atom_name == "O"))
	       || ((residue_name == "ARG") && (atom_name == "O"))
	       || ((residue_name == "ASN") && (atom_name == "O"))
	       || ((residue_name == "ASP") && ((atom_name == "O") || (atom_name == "OD1") || (atom_name == "OD2")))
	       || ((residue_name == "CYS") && ((atom_name == "O") || (atom_name == "SG")))
	       || ((residue_name == "GLN") && (atom_name == "O"))
	       || ((residue_name == "GLU") && ((atom_name == "O") || (atom_name == "OE1") || (atom_name == "OE2")))
	       || ((residue_name == "GLY") && (atom_name == "O"))
	       || ((residue_name == "HIS") && ((atom_name == "O") || (atom_name == "ND1") || (atom_name == "NE2")))
	       || ((residue_name == "ILE") && (atom_name == "O"))
	       || ((residue_name == "LEU") && (atom_name == "O"))
	       || ((residue_name == "LYS") && (atom_name == "O"))
	       || ((residue_name == "MET") && ((atom_name == "O") || (atom_name == "SD")))
	       || ((residue_name == "PHE") && (atom_name == "O"))
	       || ((residue_name == "PRO") && (atom_name == "O"))
	       || ((residue_name == "SER") && ((atom_name == "O") || (atom_name == "OG")))
	       || ((residue_name == "THR") && ((atom_name == "O") || (atom_name == "OG1")))
	       || ((residue_name == "TRP") && (atom_name == "O"))
	       || ((residue_name == "TYR") && (atom_name == "OH"))
	       || ((residue_name == "VAL") && (atom_name == "O")))
	  { 
	    acceptor_list_.push_back(atom_ptr);
	  }
      }

    return Processor::CONTINUE;
  }
} // namespace BALL
