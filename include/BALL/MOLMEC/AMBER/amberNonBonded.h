// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// $Id: amberNonBonded.h,v 1.29.4.1 2005/07/29 12:37:45 amoll Exp $
//

// Molecular Mechanics: Amber force field, non-bonded component

#ifndef BALL_MOLMEC_AMBER_NONBONDED_H
#define BALL_MOLMEC_AMBER_NONBONDED_H

#ifndef BALL_COMMON_H
#	include <BALL/common.h>
#endif

#ifndef BALL_MOLMEC_PARAMETER_LENNARDJONES_H
#	include <BALL/MOLMEC/PARAMETER/lennardJones.h>
#endif

#ifndef BALL_MOLMEC_PARAMETER_POTENTIAL1210_H
#	include <BALL/MOLMEC/PARAMETER/potential1210.h>
#endif

#ifndef BALL_MOLMEC_COMMON_FORCEFIELDCOMPONENT_H
#	include <BALL/MOLMEC/COMMON/forceFieldComponent.h>
#endif

#ifndef BALL_MOLMEC_COMMON_SUPPORT_H
#	include <BALL/MOLMEC/COMMON/support.h>
#endif

namespace BALL 
{
	/**	Amber NonBonded (VdW + Electrostatic) component
			
    	\ingroup  AMBER
	*/
	class BALL_EXPORT AmberNonBonded 
		: public ForceFieldComponent
	{
		public:

		/**	@name	Constructors and Destructors	
		*/
		//@{

		BALL_CREATE(AmberNonBonded)

		/**	Default constructor.
		*/
		AmberNonBonded()
			throw();

		/**	Constructor.
		*/
		AmberNonBonded(ForceField& force_field)
			throw();

		/**	Copy constructor
		*/
		AmberNonBonded(const AmberNonBonded& amber_non_bonded)
			throw();

		/**	Destructor.
		*/
		virtual ~AmberNonBonded()
			throw();
		//@}

		/** Assignment
		*/
		//@{
		
		/** Assignment operator
		*/
		const AmberNonBonded& operator = (const AmberNonBonded& anb)
			throw();

		/** Clear method
		*/
		virtual void clear()
			throw();

		//@}
		/** Predicates
		*/
		//@{
			
		/** Equality operator
		*/
		bool operator == (const AmberNonBonded& anb)
			throw();

		//@}
		/**	@name	Setup Methods	
		*/
		//@{

		/**	Setup method.
		*/
		virtual bool setup()
			throw(Exception::TooManyErrors);

		//@}
		/**	@name	Accessors	
		*/
		//@{

		/**	Calculates and returns the component's energy.
		*/
		virtual double updateEnergy()
			throw();

		/**	Calculates and returns the component's forces.
		*/
		virtual void updateForces()
			throw();

		/**	Update the pair list.
				This method is called by the force field whenever
				 \link ForceField::update ForceField::update \endlink  is called. It is used
				to recalculate the nonbonded pair list.
		*/
		virtual void update()
			throw(Exception::TooManyErrors);

		/**	Return the electrostatic energy.
		*/
		virtual double getElectrostaticEnergy() const
			throw();

		/**	Return the Van-der-Waals energy.
		*/
		virtual double getVdwEnergy() const
			throw();

		//@}
		/**	@name Neighbourhood and Parameter calculations
		*/
		//@{

		/**	Computes the most efficient way to calculate the non-bonded atom pairs
		*/
		virtual MolmecSupport::PairListAlgorithmType
			determineMethodOfAtomPairGeneration()
			throw();

		/**	Build a vector of non-bonded atom pairs with the vdw parameters
		*/
		virtual void buildVectorOfNonBondedAtomPairs
			(const std::vector<std::pair<Atom*, Atom*> >& atom_vector,
			 const LennardJones& lennard_jones,
			 const Potential1210& hydrogen_bond)
			throw(Exception::TooManyErrors);

		//@}

		protected:

		/*_	@name	Protected Attributes	
		*/
		//_@{

		/*_	Value of the electrostatic energy
		*/
		double	electrostatic_energy_;

		/*_	Value of the vdw energy
		*/
		double	vdw_energy_;

		//_@}

		private:

		/*_	@name	Private Attributes	
		*/
		//_@{

		/*_	Vector array with all atom pairs whose distance is smaller than cut_off
		*/
		vector<LennardJones::Data>	non_bonded_;

    /*_ Vector of flags deciding whether the pair forms a hydrogen bond or a
				standard VdW interaction.
    */
    vector<char> is_hydrogen_bond_;
 
		/*_	Number of 1-4 interactions in the vector non_bonded
		*/
		Size	number_of_1_4_;	

		/*_	Number of hydrogen bond interactions in the vector non_bonded
		*/
		Size	number_of_h_bonds_;	

		/*_	Cutoff distance for non-bonded interactions
		*/
		double	cut_off_;

		/*_	Cutoff distance for vdw interactions
		*/
		double	cut_off_vdw_;

		/*_	Cuton distance for vdw interactions
		*/
		double	cut_on_vdw_;

		/*_	Cutoff distance for electrostatic interactions
		*/
		double	cut_off_electrostatic_;

		/*_	Cuton distance for electrostatic interactions
		*/
		double	cut_on_electrostatic_;

		/*_	Inverse cube of the difference of squares of cuton and cutoff for vdW.
				This value is required for the switching function
		*/
		double inverse_distance_off_on_vdw_3_;
		
		/*_	Inverse cube of the difference of squares of cuton and cutoff for eletrostatic.
				This value is required for the switching function
		*/
		double inverse_distance_off_on_electrostatic_3_;
		
		/*_	Scaling factor for vdw_1_4_interactions
		*/
		double	scaling_vdw_1_4_;

		/*_	Scaling factor for electrostatic_1_4_interactions
		*/
		double	scaling_electrostatic_1_4_;

		/*_ Flag for using constant or distance dependent dielectric constant.
        True = distance dependent
    */
    bool    use_dist_depend_dielectric_; 

    /*_ The most efficient algorithm to calculate the non-bonded atom pairs.
        {\tt BRUTE\_FORCE}: brute force: all against all\\
        {\tt HASH\_GRID}: box grid
    */
    MolmecSupport::PairListAlgorithmType  algorithm_type_;
 		
		LennardJones	van_der_waals_;

		Potential1210 hydrogen_bond_;

		//_@}

	};
} // namespace BALL

#endif // BALL_MOLMEC_AMBER_AMBERVDW_H
