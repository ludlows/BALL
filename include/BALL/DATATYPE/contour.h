// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// $Id: contour.h,v 1.17.6.3 2005/08/12 00:32:45 amoll Exp $
//

#ifndef BALL_DATATYPE_CONTOUR_H
#define BALL_DATATYPE_CONTOUR_H

#ifndef BALL_DATATYPE_REGULARDATA2D_H
#	include <BALL/DATATYPE/regularData2D.h>
#endif

#ifndef BALL_DATATYPE_CONTOURLINE_H
#	include <BALL/DATATYPE/contourLine.h>
#endif

#include <vector>


namespace BALL
{
	/** This class is intended to create and store ContourLines belonging to the same data-set.
    	\ingroup  DatatypeMiscellaneous
	*/
	template <typename T>
	class TContour
	{
		public:

		// macro for create method
		BALL_CREATE(TContour<T>)

		/** @name    Constructors and Destructors
		*/
		//@{
		/// Default constructor
		TContour(Size num_lines=0, double start=0, double end=0);

		/// Copy constructor
		TContour(const TContour& contour);

		/// Destructor
		virtual ~TContour();
		
		//@}
		/** @name    Accessors
		*/
		//@{

		/// Gives access to the next ContourLine. Returns false if we had already returned the last line.
		bool getNextContourLine(TContourLine<T>& contour);

		//@}
		/** @name    Assignment
		*/
		//@{
			
		/// Assignment operator
		TContour& operator = (const TContour& rhs);

		/**		
		*/
		void apply(TRegularData2D<T>& data);

		/// Clear method
		virtual void clear();

		/// Reset the ContourLine-counter.
		void resetCounter();

		//@}
		/** @name    Predicates
		*/
		//@{
			
		/// Equality operator
		bool operator == (const TContour& contour) const;

		//@}

		protected:
			
		std::vector< TContourLine<T> > lines_;
		Size num_lines_;
		double start_;
		double end_;
		typename std::vector<TContourLine<T> >::const_iterator it_;
		Position index_;
	};

	/** Default type
    	\ingroup  DatatypeMiscellaneous
	 */
	typedef TContour<float> Contour;

	template <typename T>
	TContour<T>::TContour(Size num_lines, double start, double end) : lines_(num_lines), num_lines_(num_lines), start_(start), end_(end), index_(0)
	{
	}

	template <typename T>
	TContour<T>::TContour(const TContour& copyTContour) : lines_(copyTContour.lines_), num_lines_(copyTContour.num_lines_), start_(copyTContour.start_), end_(copyTContour.end_), index_(copyTContour.index_)
	{
	}

	template <typename T>
	TContour<T>::~TContour()
	{
	}

	template <typename T>
	TContour<T>& TContour<T>::operator = (const TContour& rhs)
	{
		start_     = rhs.start_;
		end_       = rhs.end_;
		num_lines_ = rhs.num_lines_;
		it_        = rhs.it_;
	}

	template <typename T>
	void TContour<T>::clear()
	{
		start_     = 0;
		end_       = 0;
		num_lines_ = 0;
		lines_     = std::vector< TContourLine<T> >(0);
		index_     = 0;
	}

	template <typename T>
	bool TContour<T>::operator == (const TContour& compTContour) const
	{
		return ((start_ == compTContour.start_) && (end_ == compTContour.end_) && (lines_ == compTContour.lines_)
			&& (num_lines_ == compTContour.num_lines_) && (it_ == compTContour.it_) && (index_ == compTContour.index_));
	}

	template <typename T>
	void TContour<T>::apply(TRegularData2D<T>& data)
	{
		Position i;
		double step = (end_ - start_) / num_lines_;

		for (i=0; i<num_lines_; i++)
		{
			TContourLine<T> con(start_ + i*step);
			con.createContourLine(data);
			lines_[i]=con;
		};

		if (num_lines_ > 0)
		{
			it_ = lines_.begin();
			index_ = 0;
		};
	}

	template <typename T>
	bool TContour<T>::getNextContourLine(TContourLine<T>& cont)
	{
		if (index_<num_lines_)
			{
				cont = *it_;
				it_++;
				index_++;
				return (true);
			} else {
				return false;
			};
	}

	template <typename T>
	void TContour<T>::resetCounter()
	{
		it_ = lines_.begin();
		index_ = 0;
	}

} // namespace BALL

#endif
