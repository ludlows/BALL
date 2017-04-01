// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// $Id: atomIterator.h,v 1.20.6.1 2005/07/28 14:01:53 amoll Exp $
//

#ifndef BALL_KERNEL_ATOMITERATOR_H
#define BALL_KERNEL_ATOMITERATOR_H

#ifndef BALL_KERNEL_ATOM_H
#	include <BALL/KERNEL/atom.h>
#endif

#ifndef BALL_KERNEL_ITERATOR_H
#	include <BALL/KERNEL/iterator.h>
#endif

#ifndef BALL_KERNEL_PREDICATE_H
#	include <BALL/KERNEL/predicate.h>
#endif

namespace BALL 
{
	/** Traits class for the AtomIterators
			
    	\ingroup  KernelIterators
	*/
	class BALL_EXPORT AtomIteratorTraits
		: public CompositeIteratorTraits 
	{
		public:

		inline AtomIteratorTraits()
			throw()
			:	CompositeIteratorTraits()
		{
			predicate_ = &RTTI::getDefault<KernelPredicate<Atom> >();
		}
			
		inline AtomIteratorTraits(const AtomIteratorTraits& traits)
			throw()
			:	CompositeIteratorTraits(traits)
		{
		}
			
		AtomIteratorTraits(const Composite& composite)
			throw();
		
		AtomIteratorTraits& operator = (const AtomIteratorTraits& traits)
			throw();

		void resetPredicate()
			throw();
	};


	inline AtomIteratorTraits::AtomIteratorTraits(const Composite& composite)
		throw()
		:	CompositeIteratorTraits(composite)
	{
		predicate_ = &RTTI::getDefault<KernelPredicate<Atom> >();
	}
		
	inline AtomIteratorTraits& AtomIteratorTraits::operator = (const AtomIteratorTraits& traits)
		throw()
	{
		CompositeIteratorTraits::operator = (traits);
		return *this;
	}

	inline void AtomIteratorTraits::resetPredicate()
		throw()
	{
		predicate_ = &RTTI::getDefault<KernelPredicate<Atom> >();
	}


	/** @name Iterator type definitions for the various AtomIterators
    	\ingroup  KernelIterators
	*/
	//@{

	/** A mutable bidirectional iterator for the Atom class.
	*/
	typedef BidirectionalIterator
		<Composite, Atom, Composite::CompositeIterator, AtomIteratorTraits>
		AtomIterator;

	/** A constant bidirectional iterator for the Atom class.
	*/
	typedef ConstBidirectionalIterator
		<Composite, Atom, Composite::CompositeIterator, AtomIteratorTraits>
		AtomConstIterator;

	/// A mutable reverse bidirectional iterator for the Atom class.
	typedef std::reverse_iterator<AtomIterator> AtomReverseIterator;

	/// A constant reverse bidirectional iterator for the Atom class.
	typedef std::reverse_iterator<AtomConstIterator> AtomConstReverseIterator;
	
	//@}
} // namespace BALL

#endif // BALL_KERNEL_ATOMITERATOR_H
