// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// $Id: comparator.h,v 1.13.6.3 2005/08/11 14:53:36 amoll Exp $
//

#ifndef BALL_CONCEPT_COMPARATOR_H
#define BALL_CONCEPT_COMPARATOR_H

#ifndef BALL_COMMON_H
#	include <BALL/common.h>
#endif

namespace BALL 
{

	/**	Generic Comparator Class.
			It is used as a baseclass, providing a common interface.
			
	 	 \ingroup ConceptsMiscellaneous
	*/
	template <typename T>
	class Comparator
	{
		public:

		/**	@name Constructors and Destructors
		*/
		//@{
	
		/**	Default constructor
		*/
		Comparator()
			throw();

		/**	Destructor
		*/
		virtual ~Comparator()
			throw();
			
		//@}
		/**	@name	Predicates
		*/
		//@{

		/**	Test if two values are equal.
		*/
		virtual bool isEqual(const T& a, const T& b) const
			throw();

		/** Test if two values are unequal.
		*/
		bool isNotEqual(const T& a, const T& b) const
			throw();

		/** Test if a is less then b.
		*/
		virtual bool isLess(const T& a, const T& b) const
			throw();

		/** Test if a is less or equal.
		*/
		bool isLessOrEqual(const T& a, const T& b) const
			throw();

		/** Test if a is greater or equal.
		*/
		bool isGreaterOrEqual(const T& a, const T& b) const
			throw();

		/** Test if a is greater then b.
		*/
		bool isGreater(const T& a, const T& b) const
			throw();

		/** Compare two values.
				-1 is returned if a  < b.  \par
				0  is returned if a == b.  \par
				1  is returned if a  > b.
		*/
		int operator () (const T& a, const T& b) const
			throw();
			
		//@}
	};

	template <typename T>
	inline Comparator<T>::Comparator()
		throw()
	{
	}

	template <typename T>
	inline Comparator<T>::~Comparator()
		throw()
	{
	}

	template <class T>
	inline bool Comparator<T>::isEqual(const T& a, const T& b) const
		throw()
	{
		return (a == b);
	}

	template <class T>
	inline bool Comparator<T>::isNotEqual(const T& a, const T& b) const
		throw()
	{
		return !isEqual(a, b);
	}

	template <class T>
	inline bool Comparator<T>::isLess(const T& a, const T& b) const
		throw()
	{
		return (a < b);
	}

	template <class T>
	inline bool Comparator<T>::isLessOrEqual(const T& a, const T& b) const
		throw()
	{
		return !isLess(b, a);
	}

	template <class T>
	inline bool Comparator<T>::isGreaterOrEqual(const T& a, const T& b) const
		throw()
	{
		return !isLess(a, b);
	}

	template <class T>
	inline bool Comparator<T>::isGreater(const T& a, const T& b) const
		throw()
	{
		return isLess(b, a);
	}

	template <class T>
	inline int Comparator<T>::operator () (const T& a,const T& b) const
		throw()
	{
		if (isEqual(a, b) == true)
		{ 
			return 0;
		} 

		if (isLess(a, b) == true)
		{ 
			return -1;
		} 
		else 
		{ 
			return 1;
		}
	}

} // namespace BALL

#endif // BALL_CONCEPT_COMPARATOR_H
