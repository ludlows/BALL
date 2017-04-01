// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// $Id: conjugateGradient.h,v 1.28.6.1 2005/07/29 12:37:47 amoll Exp $ 
//

#ifndef BALL_MOLMEC_MINIMIZATION_CONJUGATEGRADIENT_H 
#define BALL_MOLMEC_MINIMIZATION_CONJUGATEGRADIENT_H 

#ifndef	BALL_MOLMEC_MINIMIZATION_ENERGYMINIMIZER_H
#	include <BALL/MOLMEC/MINIMIZATION/energyMinimizer.h>
#endif

namespace BALL 
{ 
	/**	Conjugate gradient (CG) minimizer. 
			A conjugate gradient minimizer for geometry optimisation.
			\ingroup  MolmecEnergyMinimizer
	*/
  class BALL_EXPORT ConjugateGradientMinimizer 
		: public EnergyMinimizer
	{
    public:

		/**	@name	Enums	
		*/
		//@{

		/**	The different conjugate gradient methods implemented.
				@see	updateDirection for details on the implementation and references
		*/
		enum UpdateMethod
		{
			/**	Polak-Ribiere method
			*/
			POLAK_RIBIERE,

			/**	Fletcher-Reeves method
			*/
			FLETCHER_REEVES,

			/**	Shanno
			*/
			SHANNO
		};

		//@}	
    /**	@name	Options and Defaults
    */
    //@{

    /**	Options names
    */
    struct Option
		{
      /** The initial step length used in the line search 
      */
      static const char* STEP_LENGTH; 
    };

		/**	Defaults for all options
		*/
    struct BALL_EXPORT Default
    {
      /** The initial step length used in the line search 
      */
      static double STEP_LENGTH; 
    };


    //@}
    /**	@name	Constructors and Destructors	
    */
    //@{

    BALL_CREATE(ConjugateGradientMinimizer)

    /**	Default constructor.
    */
    ConjugateGradientMinimizer();

    /**	Constructor expecting a valid force field 
    */
    ConjugateGradientMinimizer(ForceField& force_field);

    /**	Constructor expecting a valid force field and a snapshot manager
    */
    ConjugateGradientMinimizer(ForceField& force_field,SnapShotManager *ssm);

    /**	Constructor expecting a valid force field, a snapshot manager  and options 
    */
    ConjugateGradientMinimizer
			(ForceField& force_field, SnapShotManager* ssm, const Options& options);

    /**	Constructor expecting a valid force field and options 
    */
    ConjugateGradientMinimizer(ForceField& force_field, const Options& options);

    /**	Copy constructor
    */
    ConjugateGradientMinimizer(const ConjugateGradientMinimizer& rhs);

    /**	Destructor.
    */
    virtual ~ConjugateGradientMinimizer();

    //@}
    /**	@name	Assignments 
    */
    //@{

    /**	Assignment operator
    */
    const ConjugateGradientMinimizer& operator = (const ConjugateGradientMinimizer& rhs);

    //@}
    /**	@name	Setup methods. They do all necessary preparations. 
    */
    //@{

    /**	Specific setup
    */
    virtual bool specificSetup();

    //@}
    /**	@name	Accessors 
    */
    //@{

    /** Set explicitly the value for control parameter step_length_
    */
    void setStepLength(double  value);

    /** Get the current value of option step_length_
    */
    double getStepLength() const; 

		/** Set explicitly the criterion used for updateDirection. If the
		  	method is not recognized, this function fails.
		*/
		void setUpdateMethod(UpdateMethod method);

		/** Returns the current method for updateDirection.
		 */
		UpdateMethod getUpdateMethod() const;

		/**	Calculate the next step.
				This method calculates direction and step length for the
				next optimization step. The step is not yet performed.
				This method first tries to determine a better solution 
				with a line search along the current search direction.
				If this fails, it resets the search direction to the negative 
				current gradient and performs a new line search.
				If this fails again, the step width is reduced by an order of
				magnitude and a last line search is performed along ths negative
				gradient.
				The best solution found is stored in <tt>lambda_</tt>.
				@return	bool <b>true</b> if the line search found an acceptable solution
				@see	EnergyMinimizer::findStep
		*/
		virtual double findStep();
			
		/**	Update the search direction.
				This method updates the search direction.
				It uses the Fletcher-Reeves conjugate gradient.
		*/
		virtual void updateDirection();

		/** Update the step size.
	
		*/
		virtual void updateStepSize(double lambda);

    /**	Minimize the energy of the system.
				This method executes at most <tt>iterations</tt> minimization steps.
				If the number of iterations is not given, the number specified in the
				options is taken.
				@param	iterations the maximum number of iterations
				@param	resume <b>true</b> to resume a previous run
				@see		EnergyMinimizer::minimize
    */
    virtual bool minimize(Size iterations = 0, bool resume = false); 


    protected:

    //@}
    /**	@name	Protected Attributes 
		*/
    //@{

    /*_ The step length used in the line search 
    */
    double step_; 

    /*_ The optimal lambda in the last line search 
    */
    double lambda_; 

    /*_ The number of movable atoms.
    */
		Size number_of_atoms_; 

		/*_	The update method used for the CG
		*/
		UpdateMethod method_;

		// Some variables that are needed for the Shanno criterion
		bool first_call_;
		std::vector<Vector3> a_i_, b_i_, p_t_, y_t_, p_i_, y_i_;
		double D_1_, D_4_;
		Size restart_frequency_;

    //@}

	};
} // end of namespace BALL

#endif // BALL_MOLMEC_MINIMIZATION_CONJUGATEGRADIENT_H 
