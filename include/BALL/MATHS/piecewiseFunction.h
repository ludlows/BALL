// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// $Id: piecewiseFunction.h,v 1.22.4.1 2005/07/29 12:37:39 amoll Exp $
//

#ifndef BALL_MATHS_PIECEWISEFUNCTION_H
#define BALL_MATHS_PIECEWISEFUNCTION_H

#ifndef BALL_COMMON_LIMITS_H
# include <BALL/COMMON/limits.h>
#endif

#ifndef BALL_COMMON_H
# include <BALL/common.h>
#endif

namespace BALL
{
	/** Every piece of a piecewise function needs a bunch of coefficients 
  \ingroup FunctionClasses
	*/
	typedef std::vector<double> Coefficients;

	/** An interval is defined by its limits 
	*/
	typedef std::pair<double,double> Interval;
	#ifdef INFINITY
	#undef INFINITY
	#endif
	static const double INFINITY = Limits<double>::max();

	/** Piecewise function object.
			This class provides the interface for piecewise functions needed as
			representation of radial distribution functions (@see RadialDistributionFunction). 
			It implements the <tt>Function</tt> interface. \par
			Note that intervals <b>  must </b> be disjunct and interval limits have to
			meet. We require the intervals to be sorted such that the lowest interval
			limit is the first interval of the vector. \par
			
	 */
	class BALL_EXPORT PiecewiseFunction
	{
		public:

		BALL_CREATE(PiecewiseFunction)

		/** @name Constructors and Destructors 
		*/
		//@{

		/** Default constructor 
		*/
		PiecewiseFunction() throw();

		/** Copy constructor 
		*/
		PiecewiseFunction(const PiecewiseFunction& function) throw();

		/** Detailed constructor 
		*/
		PiecewiseFunction(const std::vector<Interval>& intervals, 
			const std::vector<Coefficients>& coeffs) throw();

		/** Destructor 
		*/
		virtual ~PiecewiseFunction() throw();

		//@}
		/** @name Assignment 
		*/
		//@{

		/** Assignemnt operator 
		*/
		PiecewiseFunction& operator = (const PiecewiseFunction& function) throw();

		/** Clear function 
		*/
		void clear() throw();

		//@}
		/** @name Accessors 
		*/
		//@{

		/** Set the intervals for the piecewise definition.
				Note that this method does <b>not</b> check the definition of the
				intervals for sanity.
		 */
		void setIntervals(const std::vector<Interval>& intervals) throw();

		/** Get all the intervals 
		*/
		const std::vector<Interval>& getIntervals() const throw();

		/** Get the interval a given x belongs to 
		*/
		const Interval& getInterval(double x) const
			throw(Exception::OutOfRange);

		/** Get interval limits by index 
		*/
		const Interval& getInterval(Position index) const 
			throw(Exception::IndexOverflow);

		/** Get the interval index for a given x 
		*/
		Position getIntervalIndex(double x) const 
			throw(Exception::OutOfRange);

		/** Return the range of the definition 
		*/
		const Interval& getRange() const throw();

		/** Set the coefficients.
				Note that this method does <b>not</b> check the vector of coefficients
				for sanity.
		 */
		void setCoefficients(const vector<Coefficients>& coefficients) throw();

		/** */
		const std::vector<Coefficients>& getCoefficients() const throw();

		/** Get the coefficients for a given x 
		*/
		const Coefficients& getCoefficients(double x) const 
			throw(Exception::OutOfRange);

		/** get coefficients from index 
		*/
		const Coefficients& getCoefficients(Position index) const 
			throw(Exception::IndexOverflow);
		
		/** compute the value of the piecewise function data for a given x 
		*/
		virtual double operator () (double x) const throw();

		/** */
		void set(const std::vector<Interval>& intervals,
				const std::vector<Coefficients>& coeffs) throw();

		//@}
		/** @name Predicates 
		*/
		//@{

		/** Check whether a given x is in the range of definition 
		*/
		bool isInRange(double x) const throw();

		/** check validity of the definition 
		*/
		virtual bool isValid() const throw();

		/** Equality operator 
		*/
		bool operator == (const PiecewiseFunction& function) const throw();

		//@}
		/** @name Debugging and Diagnostics 
		*/
		//@{

		/** 	Dumps the whole content of the object 
		*/
		virtual void dump (std::ostream& s = std::cout, Size depth = 0) const throw();

		//@}

		protected:

		/*_ This vector contains the intervals of the representation 
		*/
		std::vector<Interval> intervals_;

		/*_ This vector stores the coefficients for each interval 
		*/
		std::vector<Coefficients> coefficients_;

		bool valid_;


		private:

		/*_ The range of the defnition, needed for isInRange() and getRange() 
		*/
		Interval range_;

		/*_ Set the internal range fields 
		*/
		void calculateRange() throw();

	};
}

#endif
