// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// $Id: peak.h,v 1.18.6.2 2005/08/12 14:59:09 amoll Exp $
//

#ifndef BALL_NMR_PEAK_H
#define BALL_NMR_PEAK_H

#ifndef BALL_MATHS_VECTOR3_H
#	include	<BALL/MATHS/vector3.h>
#endif

#ifndef BALL_MATHS_VECTOR2_H
#	include	<BALL/MATHS/vector2.h>
#endif

#ifndef BALL_CONCEPT_PROPERTY_H
#	include	<BALL/CONCEPT/property.h>
#endif

#include <iostream>

namespace BALL 
{
	class Atom;

	/**	Generic Peak Class.
			Each peak contains a pointer to an associated atom 
			(in the case of NMR: the atom that causes this peak).	
			 \par
			\ingroup Spectra		
	*/
	template <typename PositionType>
	class Peak
		:	public PropertyManager
	{
		public:

		/**	@name	Typedefs
		*/
		//@{
		// Type describing the coordinates and width of the peak in all its dimensions
		typedef PositionType Position;
		//@}

		/** @name	Constructors and Destructors
		*/
		//@{

		/**	Default Constructor
		*/
		Peak();
		
		/**	Copy Constructor
		*/
		Peak(const Peak& peak);
		
		/**	Destructor
		*/
		virtual ~Peak()
			throw();
		
		//@}
		/** @name Accessors
		*/
		//@{

		/** Return the peak position.
		*/
		const Position& getPosition() const;

		/** Return the peak width.
		*/
		const Position& getWidth() const;
		
		/** Return the peak intensity (amplitude).
		*/
		float getIntensity() const;
		
		/** Set the peak position.
		*/
		void setPosition(const Position& position);

		/** Set the peak width
		*/
		void setWidth(const Position& width);
		
		/** Set the peak height
		*/
		void setIntensity(float intensity);

		/**	Return the atom pointer.
		*/
		const Atom* getAtom() const;

		/**	Set the atom pointer.
		*/
		void setAtom(const Atom* atom);

		//@}
		/**	@name Assignment
		*/
		//@{

		/** Assignment Operator
		*/
		void operator = (const Peak& peak);

		//@}
		/**	@name Predicates
		*/
		//@{

		/**	Equality operator
		*/
		bool operator == (const Peak<PositionType>& peak) const;

		/**	Lesser than operator
		*/
		bool operator < (const Peak<PositionType>& peak) const;

		/**	Greater than operator
		*/
		bool operator > (const Peak<PositionType>& peak) const;
		//@}

		protected:

		Position		position_;
		Position		width_;
		float				intensity_;
		const Atom*	atom_;
	};

	template <typename PositionType>
	Peak<PositionType>::Peak()
		:	PropertyManager(),
			position_(),
			width_(),
			intensity_(0),
			atom_(0)
	{
	}

	template <typename PositionType>
	Peak<PositionType>::~Peak()
		throw()
	{
	}

	template <typename PositionType>
	Peak<PositionType>::Peak(const Peak<PositionType>& peak)
		:	PropertyManager(peak),
			position_(peak.position_),
			width_(peak.width_),
			intensity_(peak.intensity_),
			atom_(peak.atom_)
	{
	}

	template <typename PositionType>
	BALL_INLINE
	const typename Peak<PositionType>::Position& Peak<PositionType>::getPosition() const
	{
		return position_;
	}

	template <typename PositionType>
	BALL_INLINE
	const typename Peak<PositionType>::Position& Peak<PositionType>::getWidth() const
	{
		return width_;
	}

	template <typename PositionType>
	BALL_INLINE
	void Peak<PositionType>::setPosition(const typename Peak<PositionType>::Position& position)
	{
		position_ = position;
	}

	template <typename PositionType>
	BALL_INLINE
	void Peak<PositionType>::setWidth(const typename Peak<PositionType>::Position& width)
	{
		width_ = width;
	}

	template <typename PositionType>
	BALL_INLINE
	float Peak<PositionType>::getIntensity() const
	{
		return intensity_;
	}

	template <typename PositionType>
	BALL_INLINE
	void Peak<PositionType>::setIntensity(float intensity)
	{
		intensity_ = intensity;
	}

	template <typename PositionType>
	BALL_INLINE
	const Atom* Peak<PositionType>::getAtom() const
	{
		return atom_;
	}

	template <typename PositionType>
	BALL_INLINE
	void Peak<PositionType>::setAtom(const Atom* atom)
	{
		atom_ = atom;
	}

	template <typename PositionType>
	void Peak<PositionType>::operator = (const Peak<PositionType>& peak)
	{
		position_ = peak.position_;
		width_ = peak.width_;
		intensity_ = peak.intensity_;
		atom_ = peak.atom_;
	}

	template <typename PositionType>
	bool Peak<PositionType>::operator == (const Peak<PositionType>& peak) const
	{
		return ((position_ == peak.position_)
						&& (width_ == peak.width_)
						&& (intensity_ == peak.intensity_)
						&& (atom_ == peak.atom_));
	}

	template <typename PositionType>
	bool Peak<PositionType>::operator < (const Peak<PositionType>& peak) const
	{
		return (position_ < peak.position_);
	}

	template <typename PositionType>
	bool Peak<PositionType>::operator > (const Peak<PositionType>& peak) const
	{
		return (position_ > peak.position_);
	}

	/**	Output operator
	*/
	template <typename PositionType>
	std::ostream& operator << (std::ostream& os, const Peak<PositionType>& peak)	
	{
		return (os << "[ peak @ " << peak.getPosition() 
							 << ": intensity = " << peak.getIntensity()
							 << ", width = " << peak.getWidth() 
							 << "] ");
	}

	/**	@name	Convenience typedefs
		\ingroup Spectra	
	*/
	//@{
	typedef Peak<float>		Peak1D;
	typedef Peak<Vector2>	Peak2D;
	typedef Peak<Vector3>	Peak3D;
	//@}

  
} // namespace BALL

#endif // BALL_NMR_PEAK_H
