// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// $Id: amberBend.h,v 1.19.4.1 2005/07/29 12:37:44 amoll Exp $
//

// Molecular Mechanics: Amber force field, bond stretch component

#ifndef BALL_MOLMEC_AMBER_AMBERBEND_H
#define BALL_MOLMEC_AMBER_AMBERBEND_H

#ifndef BALL_COMMON_H
#	include <BALL/common.h>
#endif

#ifndef BALL_MOLMEC_COMMON_FORCEFIELDCOMPONENT_H
#	include <BALL/MOLMEC/COMMON/forceFieldComponent.h>
#endif

#ifndef BALL_MOLMEC_PARAMETER_QUADRATICANGLEBEND_H
#	include <BALL/MOLMEC/PARAMETER/quadraticAngleBend.h>
#endif

#ifndef BALL_MOLMEC_COMMON_FORCEFIELD_H
#	include <BALL/MOLMEC/COMMON/forceField.h>
#endif

namespace BALL 
{
	/**	Amber bond stretch component
			
    	\ingroup  AMBER
	*/
	class BALL_EXPORT AmberBend 
		: public ForceFieldComponent
	{
		public:

		BALL_CREATE(AmberBend)

		/** @name	Constructors and Destructors	
		*/
		//@{ 

		/**	Default constructor.
		*/
		AmberBend();

		/**	Constructor.
		*/
		AmberBend(ForceField& force_field);

		/**	Copy constructor
		*/
		AmberBend(const AmberBend& amber_stretch);

		/**	Destructor.
		*/
		virtual ~AmberBend();

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
		virtual double updateEnergy();

		/**	Calculates and returns the component's forces.
		*/
		virtual void updateForces();

		//@}

		private:

		/*_	@name	Private Attributes	
		*/
		//_@{

		/*_	pointer to the array of bends
		*/	
		vector<QuadraticAngleBend::Data>	bend_;

		QuadraticAngleBend								bend_parameters_;

		//_@}
	 
	};
} // namespace BALL

#endif // BALL_MOLMEC_AMBER_AMBERBEND_H
