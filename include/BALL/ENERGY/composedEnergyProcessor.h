// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// $Id: composedEnergyProcessor.h,v 1.12.6.2 2005/07/28 13:51:17 amoll Exp $
//

#ifndef BALL_ENERGY_COMPOSEDENERGYPROCESSOR
#define BALL_ENERGY_COMPOSEDENERGYPROCESSOR

#ifndef BALL_COMMON_H
# include <BALL/common.h>
#endif

#ifndef BALL_DATATYPE_LIST_H
# include <BALL/DATATYPE/list.h>
#endif

#ifndef BALL_ENERGY_ENERGYPROCESSOR_H
# include <BALL/ENERGY/energyProcessor.h>
#endif

namespace BALL
{
	
	/** EnergyProcessorList
			\ingroup EnergyMiscellaneous
	*/
	typedef List<EnergyProcessor*> EnergyProcessorList;

	/** A class for composed energy calculations.
			Most energies are compositions of several energy contributions. To model
			this in a convenient way this class is intended to calculate and sum a
			given list of energy contributions.
			@see EnergyProcessor
    	
			\ingroup EnergyMiscellaneous
	 */	
	class BALL_EXPORT ComposedEnergyProcessor
		: public EnergyProcessor
	{

		public:

		BALL_CREATE(ComposedEnergyProcessor)

		/** @name Constructors and destructors 
		*/
		//@{

		/** Default constructor 
		*/
		ComposedEnergyProcessor() throw();

		/** Copy constructor 
		*/
		ComposedEnergyProcessor(const ComposedEnergyProcessor& composed_energy_proc) throw();

		/** Detailed constructor 
		*/
		ComposedEnergyProcessor(EnergyProcessorList proc_list) throw();

		/** Destructor 
		*/
		virtual ~ComposedEnergyProcessor() throw();

		//@}
		/** @name Assignment 
		*/
		//@{

		/** Assignment operator 
		*/
		const ComposedEnergyProcessor& operator = 
			(const ComposedEnergyProcessor& proc) throw();

		/** Clear method.
		 * 	Clears the calculated energy and the list of processors.
		*/
		virtual void clear() throw();
	
		//@}
		/** @name Processor functions 
		*/
		//@{
		
		/** Do all calculations and sum up the different energy contributions 
		*/
		virtual bool finish() throw();

		//@}
		/** @name Accessors 
		*/
		//@{

		/** Add a component to the list of EnergyProcessors 
		*/
		void addComponent(EnergyProcessor* proc) throw();

		/** remove a component from the list 
		*/
		void removeComponent(EnergyProcessor* proc) throw();

		/** Get number of added EnergyProcessors.
		*/
		Size getNumberOfEnergyProcessors() const throw();
		
		//@}
		/** @name Predicates 
		*/
		//@{

		/** Equality operator 
		*/
		bool operator == (const ComposedEnergyProcessor& proc) const throw();

		//@}

		protected:

		EnergyProcessorList components_;


		private:

		/*_ check the validity of this instance. If one of the processors in
		 * the list is invalid, then this instance of ComposedEnergyProcessor
		 * is invalid 
		*/
		void checkValidity() throw();

	};
} // namespace BALL
#endif //  BALL_ENERGY_COMPOSEDENERGYPROCESSOR
