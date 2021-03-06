// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// $Id: structureMapper.C,v 1.29 2005/01/04 13:27:19 oliver Exp $
//

#include <BALL/STRUCTURE/structureMapper.h>

#include <BALL/COMMON/limits.h>
#include <BALL/STRUCTURE/geometricProperties.h>
#include <BALL/KERNEL/PTE.h>
#include <BALL/DATATYPE/hashGrid.h>
#include <BALL/MATHS/quaternion.h>


#include <stack>
#include <vector>
#include <map>

using namespace std;

namespace BALL
{

	/* Default constructor */
	StructureMapper::StructureMapper()
	{
	}

	/* Constructor */
	StructureMapper::StructureMapper(AtomContainer& A, AtomContainer& B)
	{
		set(A, B);
	}

	/** Destructor */
	StructureMapper::~StructureMapper()
	{
		A_ = 0;
		B_ = 0;
	}

	/* Assign the two objects to be mapped */
	void StructureMapper::set(AtomContainer & A, AtomContainer & B)
	{
		A_ = &A;
		B_ = &B;
	}

	/* Calculate the root mean squared deviation */
	double StructureMapper::calculateRMSD()
	{
		// calculate bijection, if it is not already defined
		if (bijection_.size() == 0)
		{
			calculateDefaultBijection();
		}

		// check whether we have to transform each coordinate first
		// (only in case of calculation of the RMSD of a mapping)
		bool transform = (transformation_ != Matrix4x4(1.0, 0.0, 0.0, 0.0, 
																									 0.0, 1.0, 0.0, 0.0, 
																									 0.0, 0.0, 1.0, 0.0, 
																									 0.0, 0.0, 0.0, 1.0));

		// for each pair in the bijection array, calculate square deviation for each 
		// coordinate set

		double square_deviation = 0;
		Vector3 r;

		if (transform)
		{
			for(Size i = 0; i < bijection_.size(); i++)
			{
				r = transformation_ * bijection_[i].first->getPosition();
				square_deviation += r.getSquareDistance(bijection_[i].second->getPosition());
			}
		}
		else
		{
			for(Size i = 0; i < bijection_.size(); i++)
			{
				r = bijection_[i].first->getPosition();
				square_deviation += r.getSquareDistance(bijection_[i].second->getPosition());
			}
		}

		// calculate mean square deviation
		square_deviation /=(double) bijection_.size();

		// return RMSD
		return sqrt(square_deviation);
	}

	/* Calculate the transformation to map the first of two isomorphous
	   AtomContainer objects onto the second */
	bool StructureMapper::calculateTransformation()
	{
		// check whether both acs are defined
		if ((A_ == 0) ||(B_ == 0))
		{
			return false;
		}

		// check for same number of molecular fragments(or residues)
		if (countFragments_(*A_) != countFragments_(*B_))
		{
			return false;
		}

		return true;
	}

	StructureMapper::AtomBijection StructureMapper::calculateFragmentBijection
		(const vector <Fragment*>& A, const vector<Fragment*>& B)
	{
		AtomBijection pair_array;
		AtomPairType pair;

		Size minimum = (Size)min(A.size(), B.size());

		Fragment* fragment_A = 0;
		Fragment* fragment_B = 0;

		for(Size i = 0; i < minimum; i++)
		{
			fragment_A = A[i];
			fragment_B = B[i];

			// ?????: atom names should be checked for uniqueness

			//iterate over all atoms of A and compare names with atoms of B
			for(AtomIterator atom_iterator1 = fragment_A->beginAtom(); +atom_iterator1; ++atom_iterator1)
			{
				for(AtomIterator atom_iterator2 = fragment_B->beginAtom(); +atom_iterator2; ++atom_iterator2)
				{
					if ((*atom_iterator1).getName() ==(*atom_iterator2).getName())
					{
						pair.first = &(*atom_iterator1);
						pair.second = &(*atom_iterator2);
						pair_array.push_back(pair);
					}
				}
			}
		}

		return pair_array;
	}


	bool StructureMapper::mapFragments
		(const vector<Fragment*>& A,
		 const vector<Fragment*>& B, Matrix4x4* transformation, double upper_bound, double lower_bound)
	{
		StructureMapper::AtomBijection fragment_bijection = calculateFragmentBijection(A, B);

		Size size = (Size)fragment_bijection.size();

		// if no bijection could be found, return false
		if (size == 0)
		{
		  return false;
		}

		Matrix4x4 tmp_transformation = transformation_;
		StructureMapper::AtomBijection tmp_bijection = bijection_;
	  bijection_ = fragment_bijection;

		// calculate all triangles from the bijection
		Size i, j, k;
		double square_distance;
		double min_rmsd = Limits<double>::max();
		double rmsd;

		for (k = 0; k < size; k++)
		{
			for (j = 0; j < size; j++)
			{
				square_distance =	fragment_bijection[k].first->getPosition().getSquareDistance(fragment_bijection[j].first->getPosition());
				if ((j != k) &&(square_distance >(lower_bound * lower_bound))
						&&(square_distance <(upper_bound * upper_bound)))
				{
					for (i = 0; i < size; i++)
					{
						square_distance = fragment_bijection[k].first->getPosition().getSquareDistance
								(fragment_bijection[i].first->getPosition());

						if ((i != k) &&(i != j) && (square_distance > (lower_bound * lower_bound))
								&& (square_distance < (upper_bound * upper_bound)))
						{
							transformation_ = matchPoints(fragment_bijection[k].first->getPosition(),
																						fragment_bijection[j].first->getPosition(),
																						fragment_bijection[i].first->getPosition(),
																						fragment_bijection[k].second->getPosition(),
																						fragment_bijection[j].second->getPosition(),
																						fragment_bijection[i].second->getPosition());
							rmsd = calculateRMSD();
							if (rmsd < min_rmsd)
							{
								*transformation = transformation_;
								min_rmsd = rmsd;
							}
						}
					}
				}
			}
		}

		transformation_ = tmp_transformation;
		bijection_ = tmp_bijection;

		return true;
	}

	void StructureMapper::calculateDefaultBijection()
	{
		// Remember the names of A and their atom pointers.
		StringHashMap<Atom*> A_names;
		for (AtomIterator ai = A_->beginAtom(); +ai; ++ai)
		{
			A_names.insert(std::pair<String, Atom*>(ai->getFullName(Atom::ADD_VARIANT_EXTENSIONS_AND_ID), &*ai));
		}
		
		// Iterate over all atoms of B and try to find an 
		// atom in A identical names.
		bijection_.clear();
		for (AtomIterator ai = B_->beginAtom(); +ai; ++ai)
		{
			if (A_names.has(ai->getFullName(Atom::ADD_VARIANT_EXTENSIONS_AND_ID)))
			{
				// We found two matching atoms. Remember them.
				bijection_.push_back(AtomPairType(A_names[ai->getFullName(Atom::ADD_VARIANT_EXTENSIONS_AND_ID)], &*ai));

				// Throw away the hash map entry in order to avoid
				// 1:n mappings.
				A_names.erase(ai->getFullName(Atom::ADD_VARIANT_EXTENSIONS_AND_ID));
			}
		}

		// Check whether we could map anything.
		if (bijection_.size() == 0)	
		{
			// Next stage: try to map by atom name only.
			A_names.clear();
			for (AtomIterator ai = A_->beginAtom(); +ai; ++ai)
			{
				A_names.insert(std::pair<String, Atom*>(ai->getName(), &*ai));
			}
			bijection_.clear();
			for (AtomIterator ai = B_->beginAtom(); +ai; ++ai)
			{
				if (A_names.has(ai->getName()))
				{
					// We found two matching atoms. Remember them.
					bijection_.push_back(std::pair<Atom*, Atom*>
															(A_names[ai->getName()], &*ai));
					// Throw away the hash map entry in order to avoid
					// 1:n mappings.
					A_names.erase(ai->getName());
				}
			}			
		}

		// Check whether we could map anything.
		if (bijection_.size() == 0)	
		{
			// Last stage: map in order -- first atom of A onto
			// first atom of B and so on.
			AtomIterator ai(A_->beginAtom());
			AtomIterator bi(B_->beginAtom());
			for (; +ai && +bi; ++ai, ++bi)
			{
				bijection_.push_back(std::pair<Atom*, Atom*>
														 (&*ai, &*bi));
			}
		}
	}

	Size StructureMapper::countFragments_(const AtomContainer & ac) const
	{
		Size number_of_mol_fragments = 0;

		AtomContainerConstIterator it;

		for (it = ac.beginAtomContainer(); +it; ++it)
		{
			if (RTTI::isKindOf<Fragment>(*it))
			{
				number_of_mol_fragments++;
			}
		}

		return number_of_mol_fragments;
	}



	// constructor with the following properties: The transformation maps 
	// (1) the point(vector3) w1 onto the point v1 and  
	// (2) the point w2 onto the ray that starts in v1 and goes through v2
	// (3) the point w3 into the plane generated by v1, v2 and v3

#define EPSILON 0.00001
#define EPSILON2 0.00000001

	Matrix4x4 StructureMapper::matchPoints
		(const Vector3 & w1, const Vector3 & w2, const Vector3 & w3, 
		 const Vector3 & v1, const Vector3 & v2, const Vector3 & v3)
	{
		// initialize transformation matrix
		Matrix4x4 transformation(1, 0, 0, -w1.x, 0, 1, 0, -w1.y, 0, 0, 1, -w1.z, 0, 0, 0, 1);

		// Compute the translations that map v1 and w1 onto the origin 
		// and apply them to v2, v3 and w2, w3.
		Vector3 tw2(w2.x - w1.x, w2.y - w1.y, w2.z - w1.z);
		Vector3 tw3(w3.x - w1.x, w3.y - w1.y, w3.z - w1.z);

		Vector3 tv2(v2.x - v1.x, v2.y - v1.y, v2.z - v1.z);
		Vector3 tv3(v3.x - v1.x, v3.y - v1.y, v3.z - v1.z);

		double dist_v2_v1 = tv2.getSquareLength();
		double dist_w2_w1 = tw2.getSquareLength();
		double dist_w3_w1 = tw3.getSquareLength();
		double dist_v3_v1 = tv3.getSquareLength();

		// Try to remove nasty singularities arising if the first two
		// points in each point set are too close to each other:
		//   (a) ensure (v2 != v1) 
		if ((dist_v2_v1 < EPSILON2) && (dist_v3_v1 >= EPSILON2))
		{
			tv2.swap(tv3);
		}
		//   (b) ensure (w2 != w1) 
		if ((dist_w2_w1 < EPSILON2) && (dist_w3_w1 >= EPSILON2))
		{
			tw2.swap(tw3);
		}

		Vector3 rotation_axis;
		Quaternion rotation_quat;
		Matrix4x4 rotation;
		if ((tv2.getSquareLength() >= EPSILON2) && (tw2.getSquareLength() >= EPSILON2))
		{
			// calculate the rotation axis: orthogonal to tv2 and tw2
			tw2.normalize();
			tv2.normalize();

			rotation_axis = tw2 + tv2;

			if (rotation_axis.getSquareLength() < EPSILON)
			{
				// the two axes seem to be antiparallel -
				// invert the second vector
				rotation.setIdentity();
				rotation.m11 = -1.0;
				rotation.m22 = -1.0;
				rotation.m33 = -1.0;
			}
			else
			{
				// rotate around the rotation axis
				rotation_quat.set(rotation_axis.x, rotation_axis.y, rotation_axis.z, Constants::PI);

				// Compute the matrix4x4 form of the rotation and apply it to tv3,tw2,tw3
				rotation_quat.getRotationMatrix(rotation);
			}

			tw2 = rotation * tw2;
			tw3 = rotation * tw3;

			transformation = rotation * transformation;

			if ((tw3.getSquareLength() > EPSILON2) &&(tv3.getSquareLength() > EPSILON2))
			{
				tw3.normalize();
				tv3.normalize();

				Vector3 axis_w = tv2 % tw3;
				Vector3 axis_v = tv2 % tv3;

				if ((axis_v.getSquareLength() > EPSILON2) &&(axis_w.getSquareLength() > EPSILON2))
				{
					axis_v.normalize();
					axis_w.normalize();

					rotation_axis = axis_w % axis_v;

					if (rotation_axis.getSquareLength() < EPSILON2)
					{
						double scalar_prod = axis_w * axis_v;
						if (scalar_prod < 0.0)
						{
							rotation_quat.set(tv2.x, tv2.y, tv2.z, Constants::PI);
							rotation_quat.getRotationMatrix(rotation);
						}
						else
						{
							rotation.setIdentity();
						}
					}
					else
					{
						// Compute the quaternion form of the rotation that maps tw3 onto tv3
						double angle = acos(axis_w * axis_v);
						if (angle > EPSILON)
						{
							rotation_quat.set(rotation_axis.x, rotation_axis.y, rotation_axis.z, acos(axis_w * axis_v));

							// Compute the matrix4x4 form of the rotation
							// and add it to the transformation
							rotation_quat.getRotationMatrix(rotation);
						}
						else
						{
							// Use the identity matrix instead.
							rotation.setIdentity();
						}
					}

					transformation = rotation * transformation;
				}
			}
		}

		// apply the translation onto v1
		transformation.m14 += v1.x;
		transformation.m24 += v1.y;
		transformation.m34 += v1.z;

		// done
		return transformation;
	}


	Matrix4x4 StructureMapper::matchBackboneAtoms
		(const Residue& r1, const Residue& r2)
	{
		Size counter = 0;
		
		bool	got_p1_r1 = false;
		bool  got_p2_r1 = false;
		bool	got_p3_r1 = false;
		bool  got_p1_r2 = false;
		bool  got_p2_r2 = false;
		bool  got_p3_r2 = false;

		Matrix4x4 T;
 
		Vector3 p1_r1;  // Position of C_alpha atom of residue r1
		Vector3 p2_r1;  // Position of backbone N atom of residue r1
		Vector3 p3_r1;  // Position of backbone C atom of residue r1 
		Vector3 p1_r2;  // Position of C_alpha atom of residue r2
		Vector3 p2_r2;  // Position of backbone N atom of residue r2
		Vector3 p3_r2;  // Position of backbone C atom of residue r2 

		AtomConstIterator atom_it;

		// searching the backbone atoms of residue r1
		for (atom_it = r1.beginAtom(); +atom_it; ++atom_it)
		{
			if (got_p1_r1 == false && atom_it->getName() == "CA")
			{
				p1_r1 = atom_it->getPosition();
				got_p1_r1 = true;
				counter++;
			}
			if (got_p2_r1 == false && atom_it->getName() == "N")
			{
				p2_r1 = atom_it->getPosition();
				got_p2_r1 = true;
				counter++;
			}
			if (got_p3_r1 == false && atom_it->getName() == "C")
			{
				p3_r1 = atom_it->getPosition();
				got_p3_r1 = true;
				counter++;
			}
		} 

		// searching the backbone atoms of residue r2
		for (atom_it = r2.beginAtom(); +atom_it; ++atom_it)
		{
			if (got_p1_r2 == false && atom_it->getName() == "CA")
			{
				p1_r2 = atom_it->getPosition();
				got_p1_r2 = true;
				counter++;
			}
			if (got_p2_r2 == false && atom_it->getName() == "N")
			{
				p2_r2 =(*atom_it).getPosition();
				got_p2_r2 = true;
				counter++;
			}
			if (got_p3_r2 == false && atom_it->getName() == "C")
			{
				p3_r2 = atom_it->getPosition();
				got_p3_r2 = true;
				counter++;
			}
		} 

		// Backbone atoms are missing
		if (counter != 6)
		{
			// Error: Send error message	
			Log.error() << "StructureMapper::matchBackboneAtoms: missing backbone atoms" << endl;
		}
		else
		{
			T = matchPoints(p1_r1, p2_r1, p3_r1, p1_r2, p2_r2, p3_r2);
		}

		return T;
	}

	// map the i-th residue in the list l1 
  // on the i-th residue of the list l2 (the backbone atoms are matched) 
	Size StructureMapper::mapResiduesByBackbone
		(const list<Residue*>& l1, const list<Residue*>& l2)
	{
		Size counter = 0; // number of matched residues
		Matrix4x4 null;   // the null Matrix
		TransformationProcessor T;

		// Walk down both lists and map the residues.
		list<Residue*>::const_iterator list_it_l1 = l1.begin();
		list<Residue*>::const_iterator list_it_l2 = l2.begin();
		for( ; list_it_l1 != l1.end() && list_it_l2 != l2.end(); ++list_it_l1,++list_it_l2)
		{
			// Compute the transformation matching the backbone atoms of 
			// a residue of l1 onto the corresponding residue of l2.
			T.setTransformation(matchBackboneAtoms(**list_it_l1, **list_it_l2));

			// If a valid transformation is found, (i.e. T's transformation != null),
			// apply it to the residue.
			if (!(T.getTransformation().isEqual(null)))
			{
				(*list_it_l1)->apply(T);
				counter++;
			}   
		}
	
		// Return the number of successfully matched residues.
		return(counter);
	}

	vector<vector<Fragment*> >& StructureMapper::searchPattern
		(vector<Fragment*>& pattern,
		 AtomContainer& ac, double max_rmsd, double max_center_tolerance, double upper_bound, double lower_bound)
	{
		// determine number of fragments in the pattern
		Size no_of_frag = (Size)pattern.size();

		// calculate the distances of the centers of the pattern fragments 
		// and store them in the array dist_pattern

		vector<float>	pattern_distances(no_of_frag * no_of_frag);
		vector<Vector3> pattern_centers(no_of_frag);

		Size i, j;
		GeometricCenterProcessor geo_center;

		for(i = 0; i < no_of_frag; i++)
		{
			pattern[i]->apply(geo_center);
			pattern_centers[i] = geo_center.getCenter();
		}

		float distance;

		for(i = 0; i < no_of_frag; i++)
		{
			for(j = i; j < no_of_frag; j++)
			{
				distance = pattern_centers[i].getDistance(pattern_centers[j]);
				pattern_distances[i * no_of_frag + j] = distance;
				pattern_distances[j * no_of_frag + i] = distance;
			}
		}
		pattern_centers.clear();

		// determine the molecular fragments in ac 
		// and store them in an array

		AtomContainerIterator ac_it;
		vector<Fragment*> ac_fragments;

		for	(ac_it = ac.beginAtomContainer(); 
				 ac_it != ac.endAtomContainer(); ++ac_it)
		{
			if (RTTI::isKindOf<Fragment>(*ac_it))
			{
				ac_fragments.push_back(RTTI::castTo<Fragment>(*ac_it));
			}
		}

		// determine the number of fragments of the ac
		Size no_of_comp_frag = (Size)ac_fragments.size();

		// calculate the centers of the ac fragments
		vector<Vector3> ac_centers(no_of_comp_frag);

		for(i = 1; i < no_of_comp_frag; i++)
		{
			ac_fragments[i]->apply(geo_center);
			ac_centers[i] = geo_center.getCenter();
		}

		// calculate the distances of the centers of ac fragments

		vector < float >comp_frag_dist(no_of_comp_frag * no_of_comp_frag);

		for(i = 0; i < no_of_comp_frag; i++)
		{
			for(j = i; j < no_of_comp_frag; j++)
			{
				distance = ac_centers[i].getDistance(ac_centers[j]);
				comp_frag_dist[i * no_of_comp_frag + j] = distance;
				comp_frag_dist[j * no_of_comp_frag + i] = distance;
			}
		}

		ac_centers.clear();

		// calculate an array of arrays that contains the indices of potential matching fragments

		vector < vector < Size > >indices_CF(no_of_frag);
		vector < vector < Fragment * > >* result;
		bool ready = false;
		Size counter;

		result = new vector < vector < Fragment * > >;

		for(i = 0; i < no_of_frag && ready == false; i++)
		{
			for(j = 0, counter = 0; j < no_of_comp_frag; ++j)
			{
				if (ac_fragments[j]->getName() == pattern[i]->getName())
				{
					counter++;
					indices_CF[i].push_back(j);
				}
			}
			if (counter == 0)
			{
				ready = true;
			}
		}

		// search the pattern using the array of indices

		vector < Fragment * >potential_pattern(no_of_frag);
		vector < Size > indices_of_pot_pattern(no_of_frag);
		Matrix4x4 T;
		bool distances_fit;
		Size k;
		stack < Size > index_stack;

		i = 0;
		j = 0;

		while(!ready)
		{
			indices_of_pot_pattern[i] = indices_CF[i][j];
			distances_fit = true;

			for(k = 0; k < i && distances_fit == true; k++)
			{
				distance = pattern_distances[i * no_of_frag + k] -
				comp_frag_dist[indices_of_pot_pattern[i] * no_of_comp_frag + indices_of_pot_pattern[k]];
				if (distance < -max_center_tolerance || distance > max_center_tolerance)
				{
					distances_fit = false;
				}
			}

			if (distances_fit == true)
			{
				index_stack.push(j);
				i++;
				if (i == no_of_frag)
				{
					for(k = 0; k < no_of_frag; k++)
					{
						potential_pattern.push_back(ac_fragments[indices_of_pot_pattern[k]]);

						mapFragments(potential_pattern, pattern, &T, upper_bound, lower_bound);
						if (rmsd_ <= max_rmsd)
						{
							result->push_back(potential_pattern);
							potential_pattern.clear();
						}
						else
						{
							j = 0;
						}
					}
				}
				else
				{
					j++;
					if (j == indices_CF[i].size())
					{
						i--;
						j =(Size) index_stack.top() + 1;
						index_stack.pop();
					}
				}

				if ((i == 0) &&(j == indices_CF[0].size()))
				{
					ready = true;
				}
			}
		}

		return *result;
	}

	Matrix4x4 StructureMapper::mapProteins
		(Protein& P1, Protein& P2,
		 map<String, Size>& type_map,
		 Size& no_matched_ca, double& rmsd, 
		 double upper_bound, double lower_bound, double tolerance)
	{
		// calculate bounding box of protein P1
		BoundingBoxProcessor box_processor;

		P1.apply(box_processor);

		// insert positions of CA-atoms of P1 into a three-dimensional hashgrid 
		// and in the array ca_atoms

		Vector3 upper_bound_vector(upper_bound, upper_bound, upper_bound);

		 HashGrid3 < Position > grid_P1(box_processor.getLower() - upper_bound_vector,
		  																box_processor.getUpper() - box_processor.getLower() +
			  															(float) 2.0 * upper_bound_vector, upper_bound);

		AtomIterator atom_it;
	  vector < Vector3 > ca_atoms_P1;
	  vector < Position > index_ca_P1;
		Position no_ca_P1 = 0;

		for(atom_it = P1.beginAtom(); +atom_it; ++atom_it)
		{
			if (((*atom_it).getElement() == PTE[Element::C]) &&((*atom_it).getName().trim() == "CA"))
			{
				grid_P1.insert((*atom_it).getPosition(), no_ca_P1);
				no_ca_P1++;
				ca_atoms_P1.push_back((*atom_it).getPosition());
				index_ca_P1.push_back(type_map[(*atom_it).getFragment()->getName()]);
			}
		}

		// calculate bounding box of protein P2
		P2.apply(box_processor);

		// insert positions of CA-atoms of P2 into the hashgrid grid_P2
		HashGrid3 < Position > grid_P2(box_processor.getLower() - upper_bound_vector,
																		box_processor.getUpper() - box_processor.getLower() +
																		(float) 2.0 * upper_bound_vector, upper_bound);
		Vector3 tolerance_vector(2 * tolerance, 2 * tolerance, 2 * tolerance);
		HashGrid3 < Position > fine_grid_P2(box_processor.getLower() - tolerance_vector,
																				 box_processor.getUpper() - box_processor.getLower() + tolerance_vector,
																				 2 * tolerance);

		vector < Vector3 > ca_atoms_P2;
		vector < Position > index_ca_P2;
		Size no_ca_P2 = 0;

		for(atom_it = P2.beginAtom(); +atom_it; ++atom_it)
		{
			if (((*atom_it).getElement() == PTE[Element::C]) &&((*atom_it).getName().trim() == "CA"))
			{
				grid_P2.insert((*atom_it).getPosition(), no_ca_P2);
				fine_grid_P2.insert((*atom_it).getPosition(), no_ca_P2);

				no_ca_P2++;
				ca_atoms_P2.push_back((*atom_it).getPosition());
				index_ca_P2.push_back(type_map[(*atom_it).getFragment()->getName()]);
			}
		}

		// calculate triangles between CA-atoms of P2 whose edge length are larger than lower_bound
		// and smaller than upperbound and store them in a hashgrid with respect to their edge length

		Vector3 upper(upper_bound + 1, upper_bound + 1, upper_bound + 1);
		Vector3 lower(lower_bound - 1, lower_bound - 1, lower_bound - 1);

		HashGrid3 < TVector3 < Position > >triangles_P2(lower, upper - lower, tolerance);

		HashGrid3 < Position >::BoxIterator b_it1;
		HashGridBox3 < Position >::BoxIterator b_it2, b_it3;
		HashGridBox3 < Position >::DataIterator d_it1, d_it2, d_it3;
		TVector3 < Position > index_vector;
		Vector3 distance_vector;
		float square_upper = upper_bound * upper_bound;
		float square_lower = lower_bound * lower_bound;
		float distance1, distance2, distance3;

		for(b_it1 = grid_P2.beginBox(); +b_it1; ++b_it1)
		{
			for(d_it1 =(*b_it1).beginData(); +d_it1; ++d_it1)
			{
				for(b_it2 =(*b_it1).beginBox(); +b_it2; ++b_it2)
				{
					for(d_it2 =(*b_it2).beginData(); +d_it2; ++d_it2)
					{
						if ((*d_it2) !=(*d_it1))
						{
							distance1 = ca_atoms_P2[(*d_it1)].getSquareDistance(ca_atoms_P2[(*d_it2)]);
							if (distance1 < square_upper && distance1 > square_lower)
							{
								for(b_it3 =(*b_it1).beginBox(); +b_it3; ++b_it3)
								{
									for(d_it3 =(*b_it3).beginData(); +d_it3; ++d_it3)
									{
										if ((*d_it3) !=(*d_it1) &&(*d_it3) !=(*d_it2))
										{
											distance2 = ca_atoms_P2[(*d_it1)].getSquareDistance(ca_atoms_P2[(*d_it3)]);
											if (distance2 < square_upper && distance2 > square_lower)
											{
												distance3 = ca_atoms_P2[(*d_it2)].getSquareDistance(ca_atoms_P2[(*d_it3)]);
												if (distance3 < square_upper && distance3 > square_lower)
												{
													distance1 = sqrt(distance1);
													distance2 = sqrt(distance2);
													distance3 = sqrt(distance3);
													distance_vector.set(distance1, distance2, distance3);
													index_vector.set((*d_it1),(*d_it2),(*d_it3));
													triangles_P2.insert(distance_vector, index_vector);
												}
											}
										}
									}
								}
							}
						}
					}
				}
			}
		}
		// calculate all triangles between CA-Atoms of P1 and 
		// search similar triangles between CA-Atoms of P2 stored in triangles_P2

		HashGridBox3 < TVector3 < Position > >::BoxIterator b_it4;
		HashGridBox3 < TVector3 < Position > >::DataIterator d_it4;
		HashGridBox3 < TVector3 < Position > >*box;
		HashGridBox3 < Position > *ibox;
		HashGridBox3 < Position >::BoxIterator ibox_it;
		HashGridBox3 < Position >::DataIterator id_it;

		Matrix4x4 T, T_best;
		Vector3 v;
		bool matched;
		float square_tolerance;
		Size matched_ca;

		square_tolerance = 4 * tolerance * tolerance;
		no_matched_ca = 0;
		float squared_atom_dist, current_rmsd;

		for(b_it1 = grid_P1.beginBox(); +b_it1; ++b_it1)
		{
			for(d_it1 =(*b_it1).beginData(); +d_it1; ++d_it1)
			{
				for(b_it2 =(*b_it1).beginBox(); +b_it2; ++b_it2)
				{
					for(d_it2 =(*b_it2).beginData(); +d_it2; ++d_it2)
					{
						if ((*d_it2) !=(*d_it1))
						{
							distance1 = ca_atoms_P1[(*d_it1)].getSquareDistance(ca_atoms_P1[(*d_it2)]);
							if (distance1 < square_upper && distance1 > square_lower)
							{
								distance1 = sqrt(distance1);
								for(b_it3 =(*b_it1).beginBox(); +b_it3; ++b_it3)
								{
									for(d_it3 =(*b_it3).beginData(); +d_it3; ++d_it3)
									{
										if ((*d_it3) !=(*d_it1) &&(*d_it3) !=(*d_it2))
										{
											distance2 = ca_atoms_P1[*d_it1].getSquareDistance(ca_atoms_P1[*d_it3]);
											if (distance2 < square_upper && distance2 > square_lower)
											{
												distance2 = sqrt(distance2);
												distance3 = ca_atoms_P1[*d_it2].getSquareDistance(ca_atoms_P1[*d_it3]);
												if (distance3 < square_upper && distance3 > square_lower)
												{
													distance3 = sqrt(distance3);
													distance_vector.set(distance1, distance2, distance3);
													index_vector.set(*d_it1, *d_it2, *d_it3);
													box = triangles_P2.getBox(distance_vector);

													for(b_it4 = box->beginBox(); +b_it4; ++b_it4)
													{
														for(d_it4 =(*b_it4).beginData(); +d_it4; ++d_it4)
														{
															if (index_ca_P1[(*d_it1)] == index_ca_P2[(*d_it4).x] &&
																	index_ca_P1[(*d_it2)] == index_ca_P2[(*d_it4).y] &&
																	index_ca_P1[(*d_it3)] == index_ca_P2[(*d_it4).z])
															{
																T = matchPoints(ca_atoms_P1[(*d_it1)], ca_atoms_P1[(*d_it2)], ca_atoms_P1[(*d_it3)],
																								ca_atoms_P2[(*d_it4).x], ca_atoms_P2[(*d_it4).y],	ca_atoms_P2[(*d_it4).z]);

																matched_ca = 0;
																current_rmsd = 0;
																squared_atom_dist = 0;

																for(Size i = 0; i < no_ca_P1; i++)
																{
																	v = T * ca_atoms_P1[i];
																	ibox = fine_grid_P2.getBox(v);

																	if (ibox != 0)
																	{
																		matched = false;

																		for(ibox_it = ibox->beginBox(); +ibox_it && !matched; ++ibox_it)
																		{
																			for(id_it =(*ibox_it).beginData(); +id_it && !matched; ++id_it)
																			{
																				squared_atom_dist = v.getSquareDistance(ca_atoms_P2[(*id_it)]);
																				if (squared_atom_dist <= square_tolerance)
																				{
																					matched_ca++;
																					matched = true;
																					current_rmsd += squared_atom_dist;
																				}
																			}
																		}
																	}
																}

																if (matched_ca >= no_matched_ca)
																{
																	current_rmsd = sqrt(current_rmsd / matched_ca);
																	if (matched_ca == no_matched_ca)
																	{
																		if (current_rmsd < rmsd)
																		{
																			T_best = T;
																			rmsd = current_rmsd;
																		}
																	}
																	else
																	{
																		T_best = T;
																		rmsd = current_rmsd;
																	}
																	no_matched_ca = matched_ca;
																}
															}
														}
													}
												}
											}
										}
									}
								}
							}
						}
					}
				}
			}
		}

		return T_best;
	}

}	// namespace BALL

#undef EPSILON
#undef EPSILON2
