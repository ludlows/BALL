// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// $Id: nucleotideIterator.h,v 1.12.6.1 2005/07/28 14:01:57 amoll Exp $
//

#ifndef BALL_KERNEL_NUCLEOTIDEITERATOR_H
#define BALL_KERNEL_NUCLEOTIDEITERATOR_H

#ifndef BALL_KERNEL_ITERATOR_H
#	include <BALL/KERNEL/iterator.h>
#endif

#ifndef BALL_KERNEL_PREDICATE_H
# include <BALL/KERNEL/predicate.h>
#endif

namespace BALL 
{
	class Nucleotide;

	/** NucleotideIteratorTraits
	 	  \ingroup  KernelIterators
	*/
	class BALL_EXPORT NucleotideIteratorTraits
		: public CompositeIteratorTraits
	{
		public:

		NucleotideIteratorTraits()
			:	CompositeIteratorTraits()
		{
			predicate_ = &RTTI::getDefault<KernelPredicate<Nucleotide> >();
		}
			
		NucleotideIteratorTraits(const Composite& composite)
			:	CompositeIteratorTraits(composite)
		{
			predicate_ = &RTTI::getDefault<KernelPredicate<Nucleotide> >();
		}
			
		NucleotideIteratorTraits(const NucleotideIteratorTraits& traits, bool /* deep */ = true)
			:	CompositeIteratorTraits(traits)
		{
		}
			
		NucleotideIteratorTraits& operator =(const NucleotideIteratorTraits& traits)
		{
			CompositeIteratorTraits::operator=(traits);
			return *this;
		}

		void resetPredicate()
		{
			predicate_ = &RTTI::getDefault<KernelPredicate<Nucleotide> >();
		}
	};



	typedef BidirectionalIterator
		<Composite, Nucleotide, Composite::CompositeIterator, NucleotideIteratorTraits>
		NucleotideIterator;

	typedef ConstBidirectionalIterator
		<Composite, Nucleotide, Composite::CompositeIterator, NucleotideIteratorTraits>
		NucleotideConstIterator;

	typedef std::reverse_iterator<NucleotideIterator> NucleotideReverseIterator;

	typedef std::reverse_iterator<NucleotideConstIterator> NucleotideConstReverseIterator;
} // namespace BALL

#endif // BALL_KERNEL_NUCLEOTIDEITERATOR_H
