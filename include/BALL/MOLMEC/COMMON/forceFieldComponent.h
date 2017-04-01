// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// $Id: forceFieldComponent.h,v 1.19.4.1 2005/07/30 21:53:02 amoll Exp $
//

// Molecular Mechanics: general force field component class

#ifndef BALL_MOLMEC_FORCEFIELDCOMPONENT_H
#define BALL_MOLMEC_FORCEFIELDCOMPONENT_H

#ifndef BALL_COMMON_H
#	include <BALL/common.h>
#endif

#ifndef BALL_DATATYPE_STRING_H
#	include	<BALL/DATATYPE/string.h>
#endif

#ifndef BALL_KERNEL_ATOM_H
#	include	<BALL/KERNEL/atom.h>
#endif

namespace BALL 
{
	class ForceField;

	/**	Generic force field component class.
			
    	\ingroup  MolmecCommon
	*/
	class BALL_EXPORT ForceFieldComponent 
	{

		friend class ForceField;

		public:

		/**	@name	Constructors and Destructors 
		*/
		//@{

		public:

		BALL_CREATE(ForceFieldComponent)

		/**	Default constructor.
		*/
		ForceFieldComponent();

		/**	Constructor.
		*/
		ForceFieldComponent(ForceField&	force_field);

		/**	Copy constructor.
				Creates a new instance of the given ForceFieldComponent.
				The created object will not be assigned to an existing force field.
		*/
		ForceFieldComponent(const ForceFieldComponent&	force_field_component);

		/**	Destructor.
		*/
		virtual ~ForceFieldComponent();

		//@} 

		/**	@name Setup method
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

		/**	Returns the parent force field.
				0 is returned, if no force field is assigned.
		*/	
		ForceField*	getForceField() const;

		/**	Sets the force field.
		*/
		void	setForceField(ForceField& force_field);

		/**	Sets the component's name.
		*/
		void	setName(const String& name);

		/**	Returns the component's name 
		*/
		String	getName() const;

		//@}	
		/**	@name	Force field calculations
		*/
		//@{

		/**	Returns the energy of the force field component
		*/
		virtual double getEnergy() const;

		/**	Returns the energy of the force field component.
				The current energy for this force field component is being
				calculated and returned in units of kJ/mol. \par
				If the component isn't assigned to a ForceField object,
				zero is returned.
		*/
		virtual double	updateEnergy();

		/**	Updates the atomic forces in the force field.
				The forces created by this ForceFieldComponent are
				calculated for each atom and updated in the corresponding 
				array (forces) of the ForceField instance this component 
				is assigned to.
		*/
		virtual void updateForces();

		/**	Regular update of the pair list and the like.
				Each component that contains updatable data structures
				(like nonbonded pair lists) should implement
				this method. It is called for each component of a force field
				each time  \link ForceField::update ForceField::update \endlink  is called.
				The default implementation does exactly nothing.
		*/
		virtual void update()
			throw(Exception::TooManyErrors);

		protected:

		//@}
		/*_	@name Protected Attributes
		*/
		//_@{

		/*_	The force field this component is registered in
		*/
		ForceField*	force_field_;

		/*_	The energy of the component
		*/
		double	energy_;
			
		private:
		
		//_@}
		/*_	@name Private Attributes
		*/
		//_@{

		/*_	The force field component name
		*/
		String	name_;
		
		//_@}

	};
} // namespace BALL

#endif // BALL_MOLMEC_FORCEFIELDCOMPONENT_H
