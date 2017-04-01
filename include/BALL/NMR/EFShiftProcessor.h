// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// $Id: EFShiftProcessor.h,v 1.22.4.1 2005/07/29 12:37:51 amoll Exp $
//

#ifndef BALL_NMR_EFSHIFTPROCESSOR_H
#define BALL_NMR_EFSHIFTPROCESSOR_H

#ifndef BALL_NMR_SHIFT_MODULE_H
#	include<BALL/NMR/shiftModule.h>
#endif

#ifndef BALL_KERNEL_EXPRESSION_H
#	include<BALL/KERNEL/expression.h>
#endif

namespace BALL 
{
	class Atom;
		
	/**	Shift assignment processor implementing the electric field effect. 
	\ingroup ShiftModulesNMR		
	*/
	class BALL_EXPORT EFShiftProcessor
		:	public ShiftModule
	{
		public:

		BALL_CREATE(EFShiftProcessor)

		/**	@name Enums and Constants
		*/
		//@{

		/**	A symbolic name for the electric field contribution to the chemical shift.
				@see ShiftModule::PROPERTY__SHIFT
		*/
		static const char* PROPERTY__EF_SHIFT;
		
		//@}
		/** @name	Constructors and Destructors.
		*/
		//@{

		/**	Default constructor.
		*/
		EFShiftProcessor()
			throw();
	
		/**	Copy constructor.
		*/
		EFShiftProcessor(const EFShiftProcessor& processor)
			throw();
			
		/**	Destructor.
		*/
		virtual ~EFShiftProcessor()
			throw();	

		//@}
		/**	@name	Accessors.
		*/
		//@{

		/**	Initialization method.
				This method reads the parameter section <b>ElectricFieldEffect</b> and
				parses its contents.
				This section contains the definition of two expressions that define
				a bond (the first expression matches the atom whose shift is to be
				calculated, the second describes its bond partner).
				For each of these bonds, two parameters are given, 
				$\varepsilon_1$ and	$\varepsilon_2$.  \par
				Then, this method extracts the contents of the <b>Charges</b>
				section and thus constructs a hash map containing residue and atom names 
				the corresponding charges.
				This processor is applied to all atoms in  \link operator () operator () \endlink , 
				so expect the atom charges to change!
				@see operator ()
		*/
		virtual void init()
			throw();

		//@}
		/**	@name	Processor specific functions.
		*/
		//@{
		
		/**	Processor start method.
				This method clears the bond and effector list.
				It fails if no parameters were assigned.
				@return bool, <b>false</b> if <tt>parameters_ == 0</tt>
		*/
		virtual bool start()
			throw();


		/**	operator ().
				This method sets the charge for all atoms it encounters
				(using  \link assign_charge_processor_ assign_charge_processor_ \endlink ). 
				Charged atoms	are stored in the atom list  \link effector_list_ effector_list_ \endlink .
				All bonds are stored in  \link bond_list_ bond_list_ \endlink .
				@return Processor::CONTINUE
				@param composite an arbitrary composite. All non-atom objects are ignored.
		*/
		virtual Processor::Result operator () (Composite& composite)
			throw();

		/**	Finish method.
				This method performs the chemical shift calculation.
				It iterates over all bonds stored in  \link bond_list_ bond_list_ \endlink .
				If the two bond atoms match a pair of expressions from
				 \link first_atom_expressions_ first_atom_expressions_ \endlink  and  \link second_atom_expressions_ second_atom_expressions_ \endlink ,
				the electric field vector is calculated at the bond position using
				Coulomb's law and the charges and positions of the atoms in the 
				 \link effector_list_ effector_list_ \endlink .
				The chemical shift induced by the electric field effect is calculated as
					$\delta_{EF} = \varepsilon_1 * E_z + \varepsilon_2 * E^2 $
				where constants $\varepsilon_1$ and $\varepsilon_2$ are read
				from the parameter file (section <b>ElectricFieldEffect</b>).
				The chemical shift is stored in the \emph{first} atom
				using the named property  \link ShiftModule::PROPERTY__SHIFT ShiftModule::PROPERTY__SHIFT \endlink  
				and in the named property  \link PROPERTY__EF_SHIFT PROPERTY__EF_SHIFT \endlink .
				@return bool, <b>false</b> if <tt>parameters_ == 0</tt>
		*/
		virtual bool finish()
			throw();
			
		//@}

		protected:
	
		/*_	The list of bonds collected by {\tt operator ()}.
		*/
		std::list<Bond*>				bond_list_;

		/*_	The list of charged atoms (effectors).
		*/
		std::list<Atom*>				effector_list_;

		/*_	The expressions describing the first atom of a bond.
		*/
		std::vector<Expression>	first_atom_expressions_;

		/*_	The expressions describing the first atom of a bond.
		*/
		std::vector<Expression>	second_atom_expressions_;

		/*_	The parameter $\varepsilon_1$.
		*/
		std::vector<float>			epsilon1_;

		/*_	The parameter $\varepsilon_2$.
		*/
		std::vector<float>			epsilon2_;

		/*_	The charge assignment map.
		*/
		StringHashMap<float>		charge_map_;

		/*_	A flag indicating whether effectors in the same residues are to be considered.
				Set this flag by specifying the option {\tt exclude_residue_field = true} in 
				the ElectricFieldShift section of the parameter file.
				Default is false.
		*/
		bool										exclude_residue_field_;

		/*_	A cut off value for the electric field effect.
				Any effector that is further away than this cut off is ignored.
				The distance is read from the option {\tt cut_off} in the 
				section <TT>  ElectricFieldEffect </TT> from the parameter file.
				This member contains the squared value(!) of the distance.
		*/
		float										cut_off2_;
		
 	};
  
} // namespace BALL

#endif // BALL_NMR_EFSHIFTPROCESSOR_H
