// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// $Id: steepestDescent.C,v 1.27 2005/01/25 01:05:10 amoll Exp $
//

#include <BALL/MOLMEC/MINIMIZATION/steepestDescent.h>
#include <BALL/MOLMEC/COMMON/forceField.h>
#include <BALL/MOLMEC/COMMON/snapShotManager.h>

//#define BALL_DEBUG
#undef BALL_DEBUG

using namespace std;

namespace BALL
{

	// default constructor
	 SteepestDescentMinimizer::SteepestDescentMinimizer()
		:	EnergyMinimizer(),
			line_search_(*this)
	{
	}

	// copy constructor 
	SteepestDescentMinimizer::SteepestDescentMinimizer 
		(const SteepestDescentMinimizer & minimizer)
		: EnergyMinimizer(minimizer),
			line_search_(*this)
	{
	}

	// assignment operator
	const SteepestDescentMinimizer& SteepestDescentMinimizer::operator =
		(const SteepestDescentMinimizer& minimizer)
	{
		if (&minimizer != this)
		{
			EnergyMinimizer::operator = (minimizer);
			line_search_ = minimizer.line_search_;
		}

		return *this;
	}

	// Constructor initialized with a force field
	SteepestDescentMinimizer::SteepestDescentMinimizer(ForceField& force_field)
		:	EnergyMinimizer(),
			line_search_(*this)
	{
		valid_ = setup (force_field);

		if (!valid_)
		{
			Log.error() << "SteepestDescentMinimizer: setup failed! " << std::endl;
		}
	}

	// Constructor initialized with a force field and a snapshot  manager 
	SteepestDescentMinimizer::SteepestDescentMinimizer(ForceField& force_field, SnapShotManager* ssm)
		:	EnergyMinimizer(),
			line_search_(*this)
	{
		valid_ = setup (force_field, ssm);

		if (!valid_)
		{
			Log.error() << "SteepestDescentMinimizer: setup failed! " << std::endl;
		}
	}

	// Constructor initialized with a force field and a set of options
	SteepestDescentMinimizer::SteepestDescentMinimizer
		(ForceField& force_field, const Options& new_options)
		:	EnergyMinimizer(),
			line_search_(*this)
	{
		valid_ = setup (force_field, new_options);

		if (!valid_)
		{
			Log.error() << "SteepestDescentMinimizer: setup failed! " << endl;
		}
	}

	// Constructor initialized with a force field, a snapshot manager, and a set of options
	SteepestDescentMinimizer::SteepestDescentMinimizer 
		(ForceField& force_field, SnapShotManager* ssm,
		 const Options& new_options)
		:	EnergyMinimizer(),
			line_search_(*this)
	{
		valid_ = setup(force_field, ssm, new_options);

		if (!valid_)
		{
			Log.error() << "SteepestDescentMinimizer: setup failed! " << std::endl;
		}
	}

	// destructor
	SteepestDescentMinimizer::~SteepestDescentMinimizer()
	{
	}

	// virtual function for the specific setup of derived classes
	bool SteepestDescentMinimizer::specificSetup()
	{
		maximal_number_of_iterations_ = (Size)options.setDefaultInteger
			(SteepestDescentMinimizer::Option::MAXIMAL_NUMBER_OF_ITERATIONS, 
			 (long)SteepestDescentMinimizer::Default::MAXIMAL_NUMBER_OF_ITERATIONS);

		return true;
	}

	/*  The minimizer optimizes the energy of the system 
			using a modified line search algorithm.
	*/
	bool SteepestDescentMinimizer::minimize(Size iterations, bool resume)
	{
		aborted_ = false;
		// Check for validity of minimizer and force field
		if (!isValid() || getForceField() == 0 || !getForceField()->isValid())
		{
			Log.error() << "SteepestDescentMinimizer: minimizer is not initialized correctly!" << std::endl;
			aborted_ = true;
			return false;
		}
		
		// Make sure we have something worth moving.
		if (getForceField()->getNumberOfMovableAtoms() == 0)
		{
			return true;
		}

		// If the run is to be continued, don't reset the iteration counter.
		if (!resume)
		{
			// reset the number of iterations for a restart
			setNumberOfIterations(0);
			same_energy_counter_ = 0;
		}
		Size max_iterations = std::min(getNumberOfIterations() + iterations, getMaxNumberOfIterations());
		
		LineSearch line_search(*this);
		bool converged = false;

		// Some aliases
		AtomVector& atoms(const_cast<AtomVector&>(getForceField()->getAtoms()));

		// Initial step size:
		double initial_step = 1.0;
		step_ = initial_step;

		// Compute initial energy and gradient.
		updateEnergy();
		updateForces();

		while (!converged && (getNumberOfIterations() < max_iterations))
		{
			// store the current atom positions and the current gradient/energy
			atoms.savePositions();
			storeGradientEnergy();

			// Try to find a new step
			double lambda = findStep();

			// If the result was less than zero, the line search failed.
			if (lambda < 0)
			{
				#ifdef BALL_DEBUG
					Log << "SDM: first line search failed." << std::endl;
				#endif
				// Try another step with half the initial step size
				step_ = initial_step/2.;
				lambda = findStep();
				
				// Didn't help, we abort the minimization.
				if (lambda < 0)
				{
					break;
				}
			}
	
			#ifdef BALL_DEBUG
				Log << "SDM: step taken with lambda = " << lambda << std::endl;
			#endif

			// Find a better step size for the next step
			updateStepSize(lambda);
	
			// Check for convergence.
			converged = isConverged();
			
			// End of this iteration: increment iteration counter, update pairlists etc.
			finishIteration();

			if (Maths::isNan(force_field_->getEnergy())) 
			{
				aborted_ = true;
				return false;
			}

			if (Maths::isNan(getGradient().rbegin()->x) ||
			    Maths::isNan(getGradient().rbegin()->y) ||
			    Maths::isNan(getGradient().rbegin()->z)) 
			{
				aborted_ = true;
				return false;
			}

			if (abort_by_energy_enabled_)
			{
				if (force_field_->getEnergy() > abort_energy_) 
				{
					aborted_ = true;
					return false;
				}
			}

		}
		
		return converged;
	}

	double SteepestDescentMinimizer::findStep()
	{
		// Perform a line search with the current step size.
		// We search along the current gradient: direction = current gradient
		updateDirection();

		// if the norm of the new direction is too close to zero, we assume 
		// convergence
		if (direction_.norm <= 1e-12)
		{
			return 0;
		}

		double lambda = -1.0;
		bool result = line_search_.minimize(lambda, step_ * direction_.inv_norm);

		// If the line search was not successful, return -1, else
		// the value for the optimal lambda.
		if (result == false)
		{
			return -1.0;
		}
		else
		{
			return lambda;
		}
	}

	void SteepestDescentMinimizer::updateStepSize(double lambda)
	{
		// Use the lambda we have found to make a better guess for the step size
		if (lambda < 0.2) 
		{
			step_ *= 0.5;
		} 
		else if (lambda > 0.9) 
		{
			step_ *= 2.;
		}
		#ifdef BALL_DEBUG
			Log << "SDM: new step size: " << step_ << std::endl;
		#endif
	}

	void SteepestDescentMinimizer::updateDirection()
	{
		// We search along the negative current gradient, which happens to
		// be stored in initial_gradient
		getDirection() = getInitialGradient();		
		getDirection().negate();
	}

}	// namespace BALL
