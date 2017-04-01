// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// $Id: bidirectionalIterator.h,v 1.32 2004/02/23 15:19:56 anhi Exp $ 
//

#ifndef BALL_CONCEPT_BIDIRECTIONALITERATOR_H
#define BALL_CONCEPT_BIDIRECTIONALITERATOR_H

#ifndef BALL_CONCEPT_FORWARDITERATOR_H
#	include <BALL/CONCEPT/forwardIterator.h>
#endif

namespace BALL 
{

	/**	@name	Bidirectional Iterators
	  	\ingroup  ConceptsIterators
	*/
	//@{

	/**	Constant Bidirectional Iterator.
	*/
	template <typename Container, typename DataType, typename Position, typename Traits>
	class ConstBidirectionalIterator
		: public ConstForwardIterator<Container, DataType, Position, Traits>
	{
		public:

		/** @name Typedefs
		 */
		//@{

		///
		typedef std::bidirectional_iterator_tag iterator_category;
		// convenience typedef
		typedef ConstForwardIterator<Container, DataType, Position, Traits> Base;		
		//@}

		/** @name Constructors and destructor.
		 */
		//@{

		///
		BALL_INLINE ConstBidirectionalIterator() throw() {}
	
		///
		BALL_INLINE ConstBidirectionalIterator(const ConstBidirectionalIterator& iterator) throw()
			:	Base(iterator)
		{
		}

		///
		BALL_INLINE ~ConstBidirectionalIterator() throw() {}
		//@}

		/** @name Iterator methods
		 */
		//@{

		/// Move the iterator to the beginning of the container
		BALL_INLINE void toBegin() throw(Exception::Precondition);

		/// Check whether the iterator points to the item at the beginning of the container
		BALL_INLINE bool isBegin() const throw() { return Base::getTraits().isBegin(); }

		/// Move the iterator to the position after the last item of the container
		void toEnd() throw(Exception::Precondition);

		/// Check whether the iterator points to the position after the last item of the container.
		BALL_INLINE bool isEnd() const throw() { return Base::getTraits().isEnd(); }

		/// Move the iterator to the "reverse" beginning of the container
		void toRBegin() throw(Exception::Precondition);

		/// Test whether the iterator points to the "reverse" first element	of the container
		BALL_INLINE bool isRBegin() const throw() { return Base::getTraits().isRBegin(); }

		/// Move the iterator to the position before the first element
		void toREnd()	throw(Exception::Precondition);

		/// Test wheter the iterator points to the position before the first element
		BALL_INLINE bool isREnd() const throw() { return Base::getTraits().isREnd(); }

		/// Increment operator
		BALL_INLINE ConstBidirectionalIterator& operator ++ () throw(Exception::Precondition);

		/// Postfix increment operator
		BALL_INLINE ConstBidirectionalIterator operator ++ (int) throw(Exception::Precondition);

		/// Decrement operator
		BALL_INLINE ConstBidirectionalIterator& operator -- () throw(Exception::Precondition);

		/// Postfix decrement operator
		BALL_INLINE ConstBidirectionalIterator operator -- (int) throw(Exception::Precondition);

		/// Return an iterator pointingto the first item of the container
		static ConstBidirectionalIterator begin(const Container& container) throw(Exception::Precondition);

		/// Return an iterator pointing to the position after the last element.
		static ConstBidirectionalIterator end(const Container& container) throw(Exception::Precondition);

		/// Return an iterator pointing to the last element
		static ConstBidirectionalIterator rbegin(const Container& container) throw(Exception::Precondition);

		/// Return an iterator pointing to the positon before the first element
		static ConstBidirectionalIterator rend(const Container& container) throw(Exception::Precondition);
		//@}

		protected:

		/// Construct an iterator bound to a specific container
		BALL_INLINE ConstBidirectionalIterator(const Container& container) throw()
			:	Base(container)
		{
		}
	};
	//@}

	template <typename Container, typename DataType, typename Position, typename Traits>
	BALL_INLINE
	void ConstBidirectionalIterator<Container, DataType, Position, Traits>::toBegin()
		throw(Exception::Precondition)
	{
		BALL_PRECONDITION_EXCEPTION(!Base::getTraits().isSingular(), "cannot set unbound iterator to begin")
		Base::getTraits().toBegin();
	}

	template <typename Container, typename DataType, typename Position, typename Traits>
	BALL_INLINE
	void ConstBidirectionalIterator<Container, DataType, Position, Traits>::toEnd()
		throw(Exception::Precondition)
	{
		BALL_PRECONDITION_EXCEPTION(!Base::getTraits().isSingular(), "cannot set unbound iterator to end")
		Base::getTraits().toEnd();
	}

	template <typename Container, typename DataType, typename Position, typename Traits>
	BALL_INLINE
	void ConstBidirectionalIterator<Container, DataType, Position, Traits>::toRBegin()
		throw(Exception::Precondition)
	{
		BALL_PRECONDITION_EXCEPTION(!Base::getTraits().isSingular(), "cannot set unbound iterator to reverse begin")
		Base::getTraits().toRBegin();
	}

	template <typename Container, typename DataType, typename Position, typename Traits>
	BALL_INLINE
	void ConstBidirectionalIterator<Container, DataType, Position, Traits>::toREnd()
		throw(Exception::Precondition)
	{	
		BALL_PRECONDITION_EXCEPTION(!Base::getTraits().isSingular(), "cannot set unbound iterator to reverse end")
		Base::getTraits().toREnd();
	}

	template <typename Container, typename DataType, typename Position, typename Traits>
	BALL_INLINE
	ConstBidirectionalIterator<Container, DataType, Position, Traits>& 
		ConstBidirectionalIterator<Container, DataType, Position, Traits>::operator ++ ()
		throw(Exception::Precondition)
	{
		BALL_PRECONDITION_EXCEPTION(Base::getTraits().isValid(), "cannot increment an invalid iterator")
		Base::getTraits().forward();
		return *this;
	}

	template <typename Container, typename DataType, typename Position, typename Traits>
	BALL_INLINE
	ConstBidirectionalIterator<Container, DataType, Position, Traits> 
		ConstBidirectionalIterator<Container, DataType, Position, Traits>::operator ++ (int)
		throw(Exception::Precondition)
	{
		BALL_PRECONDITION_EXCEPTION(Base::getTraits().isValid(), "cannot increment an invalid iterator")
		ConstBidirectionalIterator iterator(*this);
		++(*this);
		return iterator;
	}

	template <typename Container, typename DataType, typename Position, typename Traits>
	BALL_INLINE
	ConstBidirectionalIterator<Container, DataType, Position, Traits>& 
		ConstBidirectionalIterator<Container, DataType, Position, Traits>::operator -- ()
		throw(Exception::Precondition)
	{
		BALL_PRECONDITION_EXCEPTION(!Base::getTraits().isSingular(), "cannot decrement unbound iterator")
		Base::getTraits().backward();
		return *this;
	}

	template <typename Container, typename DataType, typename Position, typename Traits>
	BALL_INLINE
	ConstBidirectionalIterator<Container, DataType, Position, Traits> 
		ConstBidirectionalIterator<Container, DataType, Position, Traits>::operator -- (int)
		throw(Exception::Precondition)
	{
		BALL_PRECONDITION_EXCEPTION(!Base::getTraits().isSingular(), "cannot decrement an unbound iterator")
		ConstBidirectionalIterator iterator(*this);
		--(*this);
		return iterator;
	}

	template <typename Container, typename DataType, typename Position, typename Traits>
	BALL_INLINE
	ConstBidirectionalIterator<Container, DataType, Position, Traits> 
		ConstBidirectionalIterator<Container, DataType, Position, Traits>::begin(const Container& container)
		throw(Exception::Precondition)
	{
		ConstBidirectionalIterator iterator(container);
		iterator.toBegin();
		return iterator;
	}

	template <typename Container, typename DataType, typename Position, typename Traits>
	BALL_INLINE
	ConstBidirectionalIterator<Container, DataType, Position, Traits> 
		ConstBidirectionalIterator<Container, DataType, Position, Traits>::end(const Container& container)
		throw(Exception::Precondition)
	{
		ConstBidirectionalIterator iterator(container);
		iterator.toEnd();
		return iterator;
	}

	template <typename Container, typename DataType, typename Position, typename Traits>
	BALL_INLINE
	ConstBidirectionalIterator<Container, DataType, Position, Traits> 
		ConstBidirectionalIterator<Container, DataType, Position, Traits>::rbegin(const Container& container)
		throw(Exception::Precondition)
	{
		ConstBidirectionalIterator iterator(container);
		iterator.toRBegin();
		return iterator;
	}

	template <typename Container, typename DataType, typename Position, typename Traits>
	BALL_INLINE
	ConstBidirectionalIterator<Container, DataType, Position, Traits> 
		ConstBidirectionalIterator<Container, DataType, Position, Traits>::rend(const Container& container)
		throw(Exception::Precondition)
	{
		ConstBidirectionalIterator iterator(container);
		iterator.toREnd();
		return iterator;
	}

	/// Mutable bidirectional iterator
	template <typename Container, typename DataType, typename Position, typename Traits>
	class BidirectionalIterator
		: public ConstBidirectionalIterator<Container, DataType, Position, Traits>
	{
		public:

		/** @name Typedefs
		 */
		//@{
		
		///
    typedef DataType& reference;
		///
    typedef DataType* pointer;
		// convenience typedef
		typedef ConstBidirectionalIterator<Container, DataType, Position, Traits> Base;
		//@}

		/** @name Constructors and Destructor
		 */
		//@{

		/// Default constructor
		BALL_INLINE BidirectionalIterator() throw() {}
	
		/// Copy constructor
		BALL_INLINE BidirectionalIterator(const BidirectionalIterator& iterator)
			throw()
			:	ConstBidirectionalIterator<Container, DataType, Position, Traits>(iterator)
		{
		}

		/// Destructor
		BALL_INLINE ~BidirectionalIterator() throw() {}

		//@}

		/** @name Iterator methods
		 */
		//@{

		/// Dereferentiation
		BALL_INLINE reference operator * () const throw() { return (reference)Base::getTraits().getData(); }

		/// Pointer dereferentiation
		BALL_INLINE pointer operator -> () const throw() { return (pointer)&Base::getTraits().getData(); }

		/// Increment operator
		BALL_INLINE BidirectionalIterator& operator ++ ()	throw(Exception::Precondition);

		/// Postfix increment operator
		BALL_INLINE BidirectionalIterator operator ++ (int) throw(Exception::Precondition);

		/// Decrement operator
		BALL_INLINE BidirectionalIterator& operator -- ()	throw(Exception::Precondition);

		/// Postfix decrement operator
		BALL_INLINE BidirectionalIterator operator -- (int) throw(Exception::Precondition);

		/// Return an iterator pointing to the first item of the container
		static BidirectionalIterator begin(const Container& container)
			throw(Exception::Precondition);

		/// Return an iterator pointing to the position after the last element
		static BidirectionalIterator end(const Container& container)
			throw(Exception::Precondition);

		/// Return an iterator pointing to the last element.
		static BidirectionalIterator rbegin(const Container& container)
			throw(Exception::Precondition);

		/// Return an iterator pointing to the positon before the first element
		static BidirectionalIterator rend(const Container& container)
			throw(Exception::Precondition);
		//@}

		protected:

		/// Construct an iterator bound to a specific container
		BALL_INLINE BidirectionalIterator(const Container& container)	throw();
	};


	template <typename Container, typename DataType, typename Position, typename Traits>
	BALL_INLINE
	BidirectionalIterator<Container, DataType, Position, Traits>& 
		BidirectionalIterator<Container, DataType, Position, Traits>::operator ++ ()
		throw(Exception::Precondition)
	{
		Base::operator ++ ();
		return *this;
	}

	template <typename Container, typename DataType, typename Position, typename Traits>
	BALL_INLINE
	BidirectionalIterator<Container, DataType, Position, Traits> 
		BidirectionalIterator<Container, DataType, Position, Traits>::operator ++ (int)
		throw(Exception::Precondition)
	{
		BidirectionalIterator iterator(*this);
		this->operator ++ ();
		return iterator;
	}

	template <typename Container, typename DataType, typename Position, typename Traits>
	BALL_INLINE
	BidirectionalIterator<Container, DataType, Position, Traits>& 
		BidirectionalIterator<Container, DataType, Position, Traits>::operator -- ()
		throw(Exception::Precondition)
	{
		Base::operator -- ();
		return *this;
	}

	template <typename Container, typename DataType, typename Position, typename Traits>
	BALL_INLINE
	BidirectionalIterator<Container, DataType, Position, Traits> 
		BidirectionalIterator<Container, DataType, Position, Traits>::operator -- (int)
		throw(Exception::Precondition)
	{
		BidirectionalIterator iterator(*this);
		this->operator -- ();
		return iterator;
	}

	template <typename Container, typename DataType, typename Position, typename Traits>
	BALL_INLINE
	BidirectionalIterator<Container, DataType, Position, Traits> 
		BidirectionalIterator<Container, DataType, Position, Traits>::begin(const Container& container)
		throw(Exception::Precondition)
	{
		BidirectionalIterator iterator(container);
		iterator.toBegin();
		return iterator;
	}

	template <typename Container, typename DataType, typename Position, typename Traits>
	BALL_INLINE
	BidirectionalIterator<Container, DataType, Position, Traits>
  	BidirectionalIterator<Container, DataType, Position, Traits>::end(const Container& container)
		throw(Exception::Precondition)
	{
		BidirectionalIterator iterator(container);
		iterator.toEnd();
		return iterator;
	}

	template <typename Container, typename DataType, typename Position, typename Traits>
	BALL_INLINE
	BidirectionalIterator<Container, DataType, Position, Traits> 
		BidirectionalIterator<Container, DataType, Position, Traits>::rbegin(const Container& container)
		throw(Exception::Precondition)
	{
		BidirectionalIterator iterator(container);
		iterator.toRBegin();
		return iterator;
	}

	template <typename Container, typename DataType, typename Position, typename Traits>
	BALL_INLINE
	BidirectionalIterator<Container, DataType, Position, Traits> 
		BidirectionalIterator<Container, DataType, Position, Traits>::rend(const Container& container)
		throw(Exception::Precondition)
	{
		BidirectionalIterator iterator(container);
		iterator.toREnd();
		return iterator;
	}

	template <typename Container, typename DataType, typename Position, typename Traits>
	BALL_INLINE
	BidirectionalIterator<Container, DataType, Position, Traits>::BidirectionalIterator(const Container& container)
		throw()
		:	ConstBidirectionalIterator<Container, DataType, Position, Traits>(container)
	{
	}


} // namespace BALL 

#endif // BALL_CONCEPT_BIDIRECTIONALITERATOR_H
