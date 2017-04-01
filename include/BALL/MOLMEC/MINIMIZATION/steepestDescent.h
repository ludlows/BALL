// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// $Id: steepestDescent.h,v 1.23.6.2 2005/11/16 22:22:46 amoll Exp $
// Line Search Minimizer: A special class for the line search minimization algorithm

#ifndef BALL_MOLMEC_MINIMIZATION_STEEPESTDESCENT_H
#define BALL_MOLMEC_MINIMIZATION_STEEPESTDESCENT_H

#ifndef BALL_MOLMEC_MINIMIZATION_ENERGYMINIMIZER_H
# include <BALL/MOLMEC/MINIMIZATION/energyMinimizer.h>
#endif

#ifndef BALL_MOLMEC_MINIMIZATION_LINESEARCH_H
# include <BALL/MOLMEC/MINIMIZATION/lineSearch.h>
#endif

namespace BALL 
{
	/** SteepestDescentMinimizer
      \ingroup  MolmecEnergyMinimizer
	 */
	class BALL_EXPORT SteepestDescentMinimizer 
		:	public EnergyMinimizer	
	{

		public:

		/**	@name	Constructors and Destructors	
		*/
		//@{
		
		BALL_CREATE(SteepestDescentMinimizer)

		/**	Default constructor.
		*/
		SteepestDescentMinimizer();

		/**	Constructor.
		*/
		SteepestDescentMinimizer(ForceField& force_field);

		/**	Constructor.
		*/
		SteepestDescentMinimizer(ForceField& force_field,SnapShotManager *ssm);

		/**	Constructor.
		*/
		SteepestDescentMinimizer(ForceField& force_field, const Options& options);

		/**	Constructor.
		*/
		SteepestDescentMinimizer(ForceField& force_field, SnapShotManager* ssm, const Options& options);

		/**	Copy constructor
		*/
		SteepestDescentMinimizer(const SteepestDescentMinimizer& minimizer);

		/**	Destructor.
		*/
		virtual ~SteepestDescentMinimizer();

		//@}
		/**	@name	Assignments 
		*/
		//@{

		/**	Assignment operator
		*/
		const SteepestDescentMinimizer&	operator = (const SteepestDescentMinimizer& minimizer);
		//@}

		/**	@name	Setup methods 
		*/
		//@{

		/**	Specific setup
		*/
		virtual bool specificSetup();
		//@}

		/**	Minimize the energy of the system using a greedy steepest descent.
		*/
		virtual bool minimize(Size steps = 0, bool resume = false);

		/** Find the next step using a line search.
		*/
		virtual double findStep();

		/**	Update the search direction.
				Steepest descent searches along the current gradient only.
				Therefore, updateDirection only assigns direction to the last gradient
				computed (current_gradient_).
		*/
		virtual void updateDirection();

		/** Update the step size.
		 */
		virtual void updateStepSize(double lambda);
		
		protected:

		/**	The line search minimizer.
				This member is used to perform the line search in findStep
		*/
		LineSearch line_search_;

		/**	The current step size.
				This is used in findStep as an argument to the line search.
				The step size is chosen and adjusted in the  \link minimize minimize \endlink .
		*/
		double step_;
	};
} // namespace BALL

#endif // BALL_MOLMEC_MINIMIZATION_STEEPESTDESCENT_H
