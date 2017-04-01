// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// $Id: predicate.h,v 1.21.4.3 2005/08/11 23:26:19 amoll Exp $
//

#ifndef BALL_KERNEL_PREDICATE_H
#define BALL_KERNEL_PREDICATE_H

#ifndef BALL_COMMON_RTTI_H
#	include <BALL/COMMON/rtti.h>
#endif

#ifndef BALL_CONCEPT_PREDICATE_H
#	include <BALL/CONCEPT/predicate.h>
#endif

#ifndef BALL_CONCEPT_COMPOSITE_H
#	include <BALL/CONCEPT/composite.h>
#endif

namespace BALL 
{
	/**	@name	KernelPredicates Kernel predicates
			This chapter describes a set of unary predicates acting on  \link Composite \endlink instances, that
			can be used to identify the class type of a composite at runtime. \par
			This predicates are functionally similar to the  \link RTTI:isKindOf RTTI:isKindOf \endlink  function.
			They decide whether a given composite is an instance of a certain class (or of
			a class derived thereof). For instance, the application of a  \link MoleculePredicate MoleculePredicate \endlink 
			to a composite will yield <b>true</b>, if the given composite was an instance of  \link Molecule Molecule \endlink 
			(or  \link Protein Protein \endlink , as Protein is derived from Molecule). Otherwise false is returned. \par
			Each KernelPredicate possesses three methods:

				-<b>operator (const Composite&)</b>
				-<b>operator (const Composite&) const</b>
			
			The two operators are functionally identical (two methods are needed, because the class
			is derived from TUnaryPredicate and we want to make sure both methods are overwritten).
			<b>operator()</b> returns true, if the given composite is a kind of the class corresponding 
			to the respective predicate, false otherwise. \par
			
			@see	UnaryPredicate
			@see	RTTI	
    
			\ingroup  Predicates
	*/
	//@{
	

#ifndef BALL_KERNEL_PREDICATE_TYPE
#define BALL_KERNEL_PREDICATE_TYPE
	typedef UnaryPredicate<Composite>	KernelPredicateType;
#endif

	/// Kernel predicate class
	template <class T>
	class KernelPredicate
		: public UnaryPredicate<Composite>
	{
		public:

		/// Constant virtual RTTI evaluation operator
		virtual bool operator () (const Composite& composite) const
			throw();
		virtual ~KernelPredicate() {}
	};

	template <class T>
	bool KernelPredicate<T>::operator () (const Composite& composite) const
		throw()
	{
		return RTTI::isKindOf<T>(composite);
	}

	//@}

} // namespace BALL


#endif // BALL_KERNEL_PREDICATE_H
