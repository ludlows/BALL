// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// $Id: FFT2D.h,v 1.12.4.3 2005/08/12 00:32:49 amoll Exp $
//

#ifndef BALL_MATHS_TFFT2D_H
#define BALL_MATHS_TFFT2D_H

#ifndef BALL_COMMON_EXCEPTION_H
# include <BALL/COMMON/exception.h>
#endif

#ifndef BALL_DATATYPE_REGULARDATA2D_H
# include <BALL/DATATYPE/regularData2D.h>
#endif

#ifndef BALL_MATHS_VECTOR2_H
# include <BALL/MATHS/vector2.h>
#endif

#include <math.h>
#include <complex>
#include <fftw3.h>

#include <BALL/MATHS/fftwCommon.h>


namespace BALL
{
	/** A class to perform Fast Fourier Transforms and inverse Fast Fourier Transforms
			on regularly spaced two dimensional data.  \par
			This class makes use of the freely available library <b>FFTW</b>, which can be
			found at http://www.fftw.org
			coordinate system can be handled automatically. The normaliztion is chosen
			symmetrically.
			 \par
			S/TFFT2D.h
	 \ingroup FFT
	 */

	template <typename ComplexTraits>
	class TFFT2D 
		: public TRegularData2D<std::complex<typename ComplexTraits::ComplexPrecision> >
	{
		public:
		
			typedef std::complex<typename ComplexTraits::ComplexPrecision> Complex;
			typedef TRegularData2D<std::complex<typename ComplexTraits::ComplexPrecision> >	ComplexVector;

      BALL_CREATE(TFFT2D)

      /**  @name Constructors and Destructors
       */
      //@{
 
      /// Default constructor
      TFFT2D()
        throw();

			/// Copy constructor
			TFFT2D(const TFFT2D &data)
				throw();

			/** Detailed constructor.  \par
			 		@param ldnX The binary logarithm of the number of grid points in X direction (we use the logarithm to
										 ensure that the number of points is a power of two, which is important for
										 the FFT)
					@param ldnY The binary logarithm of the number of grid points in Y direction
					@param stepPhysX The step width in X direction in physical space
					@param stepPhysY The step width in Y direction in physical space
					@param origin The origin of the coordinate system
					@param inFourierSpace Flag to decide whether the data is assumed to be in physical or fourier
																space
			 */
			 // AR: ldn is not any longer the binary logarithm but the absolute number of grid points
			TFFT2D(Size ldnX, Size ldnY, double stepPhysX=1., double stepPhysY=1., Vector2 origin=Vector2(0.,0.), bool inFourierSpace=false)
				throw();
      
			/// Destructor
			virtual ~TFFT2D()
				throw();
			
			//@}

			/** @name Assignment
			 */
			//@{

			/// Assignment operator
			const TFFT2D& operator = (const TFFT2D& TFFT2D)
				throw();
			
			/** Clear the contents.
			 */
			virtual void clear()
				throw();
			
			/** Clear the contents and reset all attributes.
			 */
			virtual void destroy()
				throw();

			//@}

			/** @name Predicates
			 */
			//@{

			/** Equality operator.
			 */
			bool operator == (const TFFT2D& TFFT2D) const
				throw();
			//@}
			
			// @name Accessors
			
			/** Perform a single fast fourier transform on the data.
			 */
			void doFFT()
				throw();

			/** Perform a single inverse fourier transform on the data.
			 */
			void doiFFT()
				throw();

			/** Translate the origin in physical space about {\em trans_origin},
					i.e. the new origin will be located at the former position {\em trans_origin}.
					If the result is out of bounds, the function does nothing and
					returns <b>  false </b>.
			 */
			bool translate(const Vector2& trans_origin)
				throw();

			/** Set the step width in physical space to {\em new_width_x, new_width_y}.
				The step width in fourier space is automatically adjusted
				accordingly. {\em new_width_x and new_width_y} must be positive, otherwise
				the function does nothing and retuns <b>  false </b>.
			 */
			bool setPhysStepWidth(double new_width_x, double new_width_y)
				throw();

			/** Returns the step width in physical space in X direction.
			 */
			double getPhysStepWidthX() const
				throw();

			/** Returns the step width in physical space in Y direction.
			 */
			double getPhysStepWidthY() const
				throw();

			/** Returns the step width in fourier space in X direction.
			 */
			double getFourierStepWidthX() const
				throw();

			/** Returns the step width in fourier space in Y direction.
			 */
			double getFourierStepWidthY() const
				throw();

			/** Returns the minimal position of the grid in physical space in X direction.
			 */
			double getPhysSpaceMinX() const
				throw();

			/** Returns the minimal position of the grid in physical space in Y direction.
			 */
			double getPhysSpaceMinY() const
				throw();

			/** Returns the maximal position of the grid in physical space in X direction.
			 */
			double getPhysSpaceMaxX() const
				throw();

			/** Returns the maximal position of the grid in physical space in Y direction.
			 */
			double getPhysSpaceMaxY() const
				throw();

			/** Returns the minimal position of the grid in fourier space in X direction.
			 */
			double getFourierSpaceMinX() const
				throw();

			/** Returns the minimal position of the grid in fourier space in Y direction.
			 */
			double getFourierSpaceMinY() const
				throw();

			/** Returns the maximal position of the grid in fourier space in X direction.
			 */
			double getFourierSpaceMaxX() const
				throw();

			/** Returns the maximal position of the grid in fourier space in Y direction.
			 */
			double getFourierSpaceMaxY() const
				throw();
				
			/** AR: Return the largest grid position for the x direction. 
			 		This method returns the maximum position allowed in the grid. As the point 
					in the origin has the indices (0, 0), this method returns the number of 
					points in X direction minus one.
			  */
			Size getMaxXIndex() const
				throw();

			/** AR: Return the largest grid position for the y direction. 
			 		This method returns the maximum position allowed in the grid. As the point 
					in the origin has the indices (0, 0), this method returns the number of 
					points in Y direction minus one.
			  */
			Size getMaxYIndex() const
				throw();
				
			/** AR: Return the number of inverse transforms that have been carried out using this class.
			 		This is an important factor for the normalization of the data.
			 */
			Size getNumberOfInverseTransforms() const
				throw();

			/** AR: Returns the grid coordinate corresponding to the position.
			 */
			Vector2 getGridCoordinates(Position position) const
				throw();

			/** Returns the data at the grid position closest to <b>  pos </b>,
				and automatically includes
				the correct phase factor and (symmetric) normalization.
			 */
			Complex getData(const Vector2& pos) const
				throw(Exception::OutOfGrid);

			/** Returns the data at point <b>pos</b>. If <b>pos</b> is not a 
			 		point on the grid, the data is linearly interpolated.
					This method automatically includes the correct phase factor
					and (symmetric) normalization.
				*/
			Complex getInterpolatedValue(const Vector2& pos) const
				throw(Exception::OutOfGrid);

			/** Sets the data point at the grid position closest to <b>  pos </b>
				to the value <b>  val </b>, and -- if called in fourier space --
				automatically includes the correct phase factor and 
				(symmetric) normalization.
			 */
			void setData(const Vector2& pos, Complex val)
				throw(Exception::OutOfGrid);

			/** Access the data at the grid position closest to <b>  pos </b>.
				This function returns the "raw" data at that position.
			 */
			Complex& operator[](const Vector2& pos)
				throw(Exception::OutOfGrid);

			/** Access the data at the grid position closest to <b>  pos </b>.
			 		This function returns the "raw" data at that position.
				*/
			const Complex& operator[](const Vector2& pos) const
				throw(Exception::OutOfGrid);
				
			/** AR: Access the (raw) data at Position pos.
			 */
			Complex& operator[](const Position& pos)
				throw(Exception::OutOfGrid)
			{
				return TRegularData2D<Complex>::operator [] (pos);
			}

			/** AR: Access the (raw) data at Position pos. Const method.
				*/
			const Complex& operator[](const Position& pos) const
				throw(Exception::OutOfGrid)
			{
				return TRegularData2D<Complex>::operator [] (pos);
			}
			
			// AR:
			void setNumberOfFFTTransforms(Size num)
			{
				numPhysToFourier_ = num;
			}
			
			// AR:
			void setNumberOfiFFTTransforms(Size num)
			{
				numFourierToPhys_ = num;
			}
			
			/** This computes the phase factor in fourier space that results
				if the origin of the coordinate system in physical space
				is not in the "lower left corner".
			 */
			Complex phase(const Vector2& pos) const
				throw();
				
			/** AR: Returns <b>true</b> if the data is considered to be in Fourier space,
			 		<b>false</b> otherwise.
			 */
			bool isInFourierSpace() const
				throw();

		protected:
			Size lengthX_, lengthY_;
			bool inFourierSpace_;
			Size numPhysToFourier_;
			Size numFourierToPhys_;
			Vector2 origin_;
			double stepPhysX_, stepPhysY_;
			double stepFourierX_, stepFourierY_;
      Vector2 minPhys_, maxPhys_;
      Vector2 minFourier_, maxFourier_;
      
      // AR: new version for FFTW3
			typename ComplexTraits::FftwPlan planForward_;
			typename ComplexTraits::FftwPlan planBackward_;

			
			// AR: to control plan calculation with new fftw3
			Size dataLength_;
			Complex *dataAdress_;
			bool planCalculated_;
			
	};
	
	/**	Default type
	*/
	typedef TFFT2D<DoubleTraits> FFT2D;
	
	// AR:
	/** Global assignment operator from TFFT2D to TRegularData2D<Complex>
	 */
	template <typename ComplexTraits>
	const TRegularData2D<typename TFFT2D<ComplexTraits>::Complex>& operator<< 
			(TRegularData2D<typename TFFT2D<ComplexTraits>::Complex>& to, const TFFT2D<ComplexTraits>& from)
		throw();
	
	/** Global assignment operator from FFT3D to TRegularData3D<float>.
	 		This operator assigns the <b>real</b> part of the complex TFFT2D-data to the
			TRegularData2D<float> to.
	 */
	template <typename ComplexTraits>
	const RegularData2D& operator << (RegularData2D& to, const TFFT2D<ComplexTraits>& from)
		throw();
	
	template <typename ComplexTraits>
	TFFT2D<ComplexTraits>::TFFT2D()
		throw()
		: TRegularData2D<Complex>(),
			dataLength_(0),
			dataAdress_(0),
			planCalculated_(false)
	{
	}
	
	template <typename ComplexTraits>
	bool TFFT2D<ComplexTraits>::operator == (const TFFT2D<ComplexTraits>& fft2D) const
		throw()
	{
		// AR: test whether data_.size() == fft2D.data_.size()
		//     instead of testing 2 lengths. Better for vector handling.
		
		if (lengthX_ == fft2D.lengthX_ &&
				lengthY_ == fft2D.lengthY_ &&
				origin_ == fft2D.origin_ &&
				stepPhysX_ == fft2D.stepPhysX_ &&
				stepPhysY_ == fft2D.stepPhysY_ &&
				stepFourierX_ == fft2D.stepFourierX_ &&
				stepFourierY_ == fft2D.stepFourierY_ &&
				minPhys_ == fft2D.minPhys_ &&
				maxPhys_ == fft2D.maxPhys_ &&
				minFourier_ == fft2D.minFourier_ &&
				maxFourier_ == fft2D.maxFourier_ &&
				numPhysToFourier_ == fft2D.numPhysToFourier_ &&
				numFourierToPhys_ == fft2D.numFourierToPhys_)
		{
			Vector2 min  = inFourierSpace_ ?  minFourier_  :   minPhys_;
			Vector2 max  = inFourierSpace_ ?  maxFourier_  :   maxPhys_;
			double stepX = inFourierSpace_ ? stepFourierX_ : stepPhysX_;
			double stepY = inFourierSpace_ ? stepFourierY_ : stepPhysY_;	
			
			for (double posX=min.x; posX<=max.x; posX+=stepX)
			{
				for (double posY=min.y; posY<=max.y; posY+=stepY)
				{
					if (getData(Vector2(posX,posY)) != fft2D.getData(Vector2(posX,posY)))
					{
						return false;
					}
				}
			}
			
			return true;
		}
	
		return false;
	}
	
	template <typename ComplexTraits>
	bool TFFT2D<ComplexTraits>::translate(const Vector2& trans_origin)
		throw()
	{
		Position internalOriginX = (Position) rint(trans_origin.x*stepPhysX_);
		Position internalOriginY = (Position) rint(trans_origin.y*stepPhysY_);
		
		if ((internalOriginX <= lengthX_) && (internalOriginY <= lengthY_))
		{
			origin_.x = trans_origin.x;
			origin_.y = trans_origin.y;
			
      minPhys_ = Vector2(-origin_.x,-origin_.y);
  	  maxPhys_ = Vector2(((lengthX_-1)*stepPhysX_)-origin_.x,((lengthY_-1)*stepPhysY_)-origin_.y);
	 	  minFourier_ = Vector2(-(lengthX_/2.-1)*stepFourierX_,-(lengthY_/2.-1)*stepFourierY_);
			maxFourier_ = Vector2((lengthX_/2.)*stepFourierX_,(lengthY_/2.)*stepFourierY_);
		 
			return true;
		}
		else
		{
			return false;
		}
	}

	template <typename ComplexTraits>
	bool TFFT2D<ComplexTraits>::setPhysStepWidth(double new_width_x, double new_width_y)
		throw()
	{
		if ((new_width_x <= 0) || (new_width_y <= 0))
		{
			return false;
		}
		else
		{
			stepPhysX_ = new_width_x;
			stepPhysY_ = new_width_y;
			stepFourierX_ = 2.*M_PI/(stepPhysX_*lengthX_);
			stepFourierY_ = 2.*M_PI/(stepPhysY_*lengthY_);

			minPhys_ = Vector2(-origin_.x,-origin_.y);
  	  maxPhys_ = Vector2(((lengthX_-1)*stepPhysX_)-origin_.x,((lengthY_-1)*stepPhysY_)-origin_.y);
	 	  minFourier_ = Vector2(-(lengthX_/2.-1)*stepFourierX_,-(lengthY_/2.-1)*stepFourierY_);
			maxFourier_ = Vector2((lengthX_/2.)*stepFourierX_,(lengthY_/2.)*stepFourierY_);
	
      return true;
		}
	}
	
	template <typename ComplexTraits>
	double TFFT2D<ComplexTraits>::getPhysStepWidthX() const
		throw()
	{
		return stepPhysX_;
	}

	template <typename ComplexTraits>
	double TFFT2D<ComplexTraits>::getPhysStepWidthY() const
		throw()
	{
		return stepPhysY_;
	}

	template <typename ComplexTraits>
	double TFFT2D<ComplexTraits>::getFourierStepWidthX() const
		throw()
	{
		return stepFourierX_;
	}

	template <typename ComplexTraits>
	double TFFT2D<ComplexTraits>::getFourierStepWidthY() const
		throw()
	{
		return stepFourierY_;
	}

	template <typename ComplexTraits>
	double TFFT2D<ComplexTraits>::getPhysSpaceMinX() const
		throw()
	{
    return minPhys_.x;
  }

	template <typename ComplexTraits>
	double TFFT2D<ComplexTraits>::getPhysSpaceMinY() const
		throw()
	{
    return minPhys_.y;
  }

	template <typename ComplexTraits>
	double TFFT2D<ComplexTraits>::getPhysSpaceMaxX() const
		throw()
	{
    return maxPhys_.x;
 	}

	template <typename ComplexTraits>
	double TFFT2D<ComplexTraits>::getPhysSpaceMaxY() const
		throw()
	{
    return maxPhys_.y;
 	}

	template <typename ComplexTraits>
	double TFFT2D<ComplexTraits>::getFourierSpaceMinX() const
		throw()
	{
		return minFourier_.x;
	}

	template <typename ComplexTraits>
	double TFFT2D<ComplexTraits>::getFourierSpaceMinY() const
		throw()
	{
		return minFourier_.y;
	}

	template <typename ComplexTraits>
	double TFFT2D<ComplexTraits>::getFourierSpaceMaxX() const
		throw()
	{
		return maxFourier_.x;
	}

	template <typename ComplexTraits>
	double TFFT2D<ComplexTraits>::getFourierSpaceMaxY() const
		throw()
	{
		return maxFourier_.y;
	}
	
	template <typename ComplexTraits>
	Size TFFT2D<ComplexTraits>::getMaxXIndex() const
		throw()
	{
		return (lengthX_ - 1);
	}
	
	template <typename ComplexTraits>
	Size TFFT2D<ComplexTraits>::getMaxYIndex() const
		throw()
	{
		return (lengthY_ - 1);
	}
	
	template <typename ComplexTraits>
	Size TFFT2D<ComplexTraits>::getNumberOfInverseTransforms() const
		throw()
	{
		return numFourierToPhys_;
	}
	
	// AR: new for compatibility with FFT3D
	template <typename ComplexTraits>
	Vector2 TFFT2D<ComplexTraits>::getGridCoordinates(Position position) const
		throw()
	{
		if (!inFourierSpace_)
		{
			if (position >= ComplexVector::size())
			{
				throw Exception::OutOfGrid(__FILE__, __LINE__);
			}
		
			Vector2 r;
			Position  x, y;


			// AR: ??????
			y = position % lengthY_;
			x = position / lengthY_;

			r.set(-origin_.x + (float)x * stepPhysX_,
						-origin_.y + (float)y * stepPhysY_);

			return r;
		}
		else
		{
			if (position >= ComplexVector::size())
			{
				throw Exception::OutOfGrid(__FILE__, __LINE__);
			}
		
			Vector2 r;
			Index x, y;
	
			// AR: ??????
			y = position % lengthY_;
			x = position / lengthY_;

			if (x>=lengthX_/2.)
			{
				x-=lengthX_;
			}
			
			if (y>=lengthY_/2.)
			{
				y-=lengthY_;
			}

			r.set((float)x * stepFourierX_,
						(float)y * stepFourierY_);

			return r;
		}
	}
	
	
	
	template <typename ComplexTraits>
	typename TFFT2D<ComplexTraits>::Complex TFFT2D<ComplexTraits>::getData(const Vector2& pos) const
		throw(Exception::OutOfGrid)
	{
		Complex result;
		double normalization=1.;

		if (!inFourierSpace_)
		{
			result = (*this)[pos];
			normalization=1./((float)pow((float)(lengthX_*lengthY_),(int)numFourierToPhys_));
		}
		else
		{
			result = (*this)[pos] * phase(pos);
			normalization=1./(2.*M_PI)*(stepPhysX_*stepPhysY_)/((float)pow((float)(lengthX_*lengthY_),(int)numFourierToPhys_));
		}

		result *= normalization;
		
		return result;
	}

	template <typename ComplexTraits>
	typename TFFT2D<ComplexTraits>::Complex TFFT2D<ComplexTraits>::getInterpolatedValue(const Vector2& pos) const
		throw(Exception::OutOfGrid)
	{
		Complex result;
		
		Vector2 min  = inFourierSpace_ ? minFourier_   :   minPhys_;
		Vector2 max  = inFourierSpace_ ? maxFourier_   :   maxPhys_;
		double stepX = inFourierSpace_ ? stepFourierX_ : stepPhysX_;
		double stepY = inFourierSpace_ ? stepFourierY_ : stepPhysY_;
		
		if (    (pos.x < min.x) || (pos.y < min.y)
				 || (pos.x > max.x) || (pos.y > max.y)  )
		{
			throw Exception::OutOfGrid(__FILE__, __LINE__);
		}

		Vector2 h(pos.x - min.x, pos.y - min.y);
		double modX = fmod((double)h.x,stepX);
		double modY = fmod((double)h.y,stepY);

		if (modX==0 && modY ==0) // we are on the grid
		{
			return getData(pos);
		}

		double beforeX = floor(h.x/stepX)*stepX+ min.x;
		double beforeY = floor(h.y/stepY)*stepY+ min.y;
		double afterX  =  ceil(h.x/stepX)*stepX+ min.x;
		double afterY  =  ceil(h.y/stepY)*stepY+ min.y;
			
		double tx = (pos.x - beforeX)/stepX;
		double ty = (pos.y - beforeY)/stepY;

		result  = getData(Vector2(beforeX,beforeY))*(typename ComplexTraits::ComplexPrecision)((1.-tx)*(1.-ty));
		result += getData(Vector2(afterX, beforeY))*(typename ComplexTraits::ComplexPrecision)(    tx *(1.-ty));
		result += getData(Vector2(beforeX,afterY ))*(typename ComplexTraits::ComplexPrecision)((1.-tx)*    ty );
		result += getData(Vector2(afterX, afterY ))*(typename ComplexTraits::ComplexPrecision)(    tx *    ty );

		return result;
	}

	template <typename ComplexTraits>
	void TFFT2D<ComplexTraits>::setData(const Vector2& pos, Complex val)
		throw(Exception::OutOfGrid)
	{
		Complex dummy;
	
		if (!inFourierSpace_)
		{
			dummy = Complex(val.real()*((float)pow((float)(lengthX_*lengthY_),(int)numFourierToPhys_)),
												val.imag()*((float)pow((float)(lengthX_*lengthY_),(int)numFourierToPhys_)));
	
			(*this)[pos]=dummy;
		}
		else
		{
			val*=phase(pos)*(typename ComplexTraits::ComplexPrecision)((2*M_PI/(stepPhysX_*stepPhysY_)))
										 *(typename ComplexTraits::ComplexPrecision)pow((typename ComplexTraits::ComplexPrecision)(lengthX_*lengthY_),(int)numFourierToPhys_);
			
			dummy = val;
		
			(*this)[pos]=dummy;
		}
	}

	template <typename ComplexTraits>
	typename TFFT2D<ComplexTraits>::Complex& TFFT2D<ComplexTraits>::operator[](const Vector2& pos)
		throw(Exception::OutOfGrid)
	{
		Index internalPos;

		if (!inFourierSpace_)
		{
			Index i, j;
			
			i = (Index) rint((pos.x+origin_.x)/stepPhysX_);
			j = (Index) rint((pos.y+origin_.y)/stepPhysY_);

			internalPos = j + i*lengthY_;
		}
		else
		{
			Index i, j;

			i = (Index) rint(pos.x/stepFourierX_);
			j = (Index) rint(pos.y/stepFourierY_);

			if (i<0)
			{
				i+=lengthX_;
			}

			if (j<0)
			{
				j+=lengthY_;
			}

			internalPos = (j + i*lengthY_);
		}

		if ((internalPos < 0) || (internalPos>=(Index) (lengthX_*lengthY_)))
		{
			throw Exception::OutOfGrid(__FILE__, __LINE__);
		}
		
		return operator [] (internalPos);
	}

	template <typename ComplexTraits>
	const typename TFFT2D<ComplexTraits>::Complex& TFFT2D<ComplexTraits>::operator[](const Vector2& pos) const
		throw(Exception::OutOfGrid)
	{
		Index internalPos;

		if (!inFourierSpace_)
		{
			Index i, j;
			
			i = (Index) rint((pos.x+origin_.x)/stepPhysX_);
			j = (Index) rint((pos.y+origin_.y)/stepPhysY_);

			internalPos = j + i*lengthY_;
		}
		else
		{
			Index i, j;

			i = (Index) rint(pos.x/stepFourierX_);
			j = (Index) rint(pos.y/stepFourierY_);

			if (i<0)
			{
				i+=lengthX_;
			}

			if (j<0)
			{
				j+=lengthY_;
			}

			internalPos = (j + i*lengthY_);
		}

		if ((internalPos < 0) || (internalPos>=(Index) (lengthX_*lengthY_)))
		{
			throw Exception::OutOfGrid(__FILE__, __LINE__);
		}
		
		return operator [] (internalPos);
	}
	
	template <typename ComplexTraits>
	typename TFFT2D<ComplexTraits>::Complex TFFT2D<ComplexTraits>::phase(const Vector2& pos) const
		throw()
	{
	  double phase = 2.*M_PI*(  rint(pos.x/stepFourierX_)*rint(origin_.x/stepPhysX_)
															/lengthX_
														+ rint(pos.y/stepFourierY_)*rint(origin_.y/stepPhysY_)
															/lengthY_ );

		Complex result = Complex(cos(phase), sin(phase));
						
		return result;
	}
	
	template <typename ComplexTraits>
	bool TFFT2D<ComplexTraits>::isInFourierSpace() const
		throw()
	{
		return inFourierSpace_;
	}
	
	template <typename ComplexTraits>
	const TRegularData2D<typename TFFT2D<ComplexTraits>::Complex>& operator << 
		(TRegularData2D<typename TFFT2D<ComplexTraits>::Complex>& to, const TFFT2D<ComplexTraits>& from)
		throw()
	{
		// first decide if the FFT3D data is in Fourier space.
		if (!from.isInFourierSpace())
		{
			// create a new grid
			Size lengthX = from.getMaxXIndex()+1;
			Size lengthY = from.getMaxYIndex()+1;
			
			TRegularData2D<typename TFFT2D<ComplexTraits>::Complex> newGrid(TRegularData2D<typename TFFT2D<ComplexTraits>::Complex>::IndexType(lengthX, lengthY),
																			Vector2(from.getPhysSpaceMinX(), from.getPhysSpaceMinY()),
																			Vector2(from.getPhysSpaceMaxX(), from.getPhysSpaceMaxY()));

			// and fill it
			double normalization=1./(pow((float)(lengthX*lengthY),from.getNumberOfInverseTransforms()));
			typename TFFT2D<ComplexTraits>::Complex dataIn;
			typename TFFT2D<ComplexTraits>::Complex dataOut;
			
			for (Position i = 0; i < from.size(); i++)
			{
				Position x, y;

				y =  i % lengthY;
				x =  i / lengthY;

				dataIn  = from[i];
				dataOut = dataIn;
				
				newGrid[x + y*lengthY] = dataOut*(typename ComplexTraits::ComplexPrecision)normalization;
			}

			to = newGrid;

			return to;
		}
		else
		{
			// we are in Fourier space
			
			// create a new grid
			Size lengthX = from.getMaxXIndex()+1;
			Size lengthY = from.getMaxYIndex()+1;
			//float stepPhysX = from.getPhysStepWidthX();
			//float stepPhysY = from.getPhysStepWidthY();
			float stepFourierX = from.getFourierStepWidthX();
			float stepFourierY = from.getFourierStepWidthY();


			
			TRegularData2D<typename TFFT2D<ComplexTraits>::Complex> newGrid(TRegularData2D<typename TFFT2D<ComplexTraits>::Complex>::IndexType(lengthX, lengthY),
																			Vector2(from.getFourierSpaceMinX(), 
																							from.getFourierSpaceMinY()),
																			Vector2(from.getFourierSpaceMaxX(),
																							from.getFourierSpaceMaxY()));

			// and fill it
			// AR: old double normalization=1./(sqrt(2.*M_PI))*(stepPhysX*stepPhysY*stepPhysZ)/(pow((float)(lengthX*lengthY*lengthZ),from.getNumberOfInverseTransforms()));
			double normalization=1./(2.*M_PI)/(pow((float)(lengthX*lengthY),from.getNumberOfInverseTransforms()));
			
			
			Index x, y;
			Vector2 r;
			typename TFFT2D<ComplexTraits>::Complex dataIn;
			typename TFFT2D<ComplexTraits>::Complex dataOut;
	
			for (Position i = 0; i < from.size(); i++)
			{
				y =  i % lengthY;
				x =  i / lengthY;

				if (x>lengthX/2.)
				{
					x-=lengthX;
				}

				if (y>lengthY/2.)
				{
					y-=lengthY;
				}

				r.set((float)x * stepFourierX,
							(float)y * stepFourierY);

				dataIn = from[i];
				dataOut = dataIn;
				
				newGrid[x + y*lengthY] = dataOut*(typename ComplexTraits::ComplexPrecision)normalization*from.phase(r);
			}

			to = newGrid;

			return to;
		}
	}
	
	template <typename ComplexTraits>
	const RegularData2D& operator << (RegularData2D& to, const TFFT2D<ComplexTraits>& from)
		throw()
	{
		// first decide if the FFT3D data is in Fourier space.
		if (!from.isInFourierSpace())
		{
			// create a new grid
			Size lengthX = from.getMaxXIndex()+1;
			Size lengthY = from.getMaxYIndex()+1;
			
			RegularData2D newGrid(RegularData2D::IndexType(lengthX, lengthY),
														Vector2(from.getPhysSpaceMinX(), 
																		from.getPhysSpaceMinY()),
														Vector2(from.getPhysSpaceMaxX(),
																		from.getPhysSpaceMaxY()));

			// and fill it
			double normalization = 1./(pow((float)(lengthX*lengthY),from.getNumberOfInverseTransforms()));
			typename TFFT2D<ComplexTraits>::Complex dataIn;
			typename TFFT2D<ComplexTraits>::Complex dataOut;
			
			for (Position i = 0; i < from.size(); i++)
			{
				Position x, y;

				y =  i % lengthY;
				x =  i / lengthY;

				dataIn  = from[i];
				dataOut = dataIn;
				
				newGrid[x + y*lengthY] = dataOut.real()*normalization;
			}

			to = newGrid;

			return to;
		}
		else
		{
			// we are in Fourier space
			
			// create a new grid
			Size lengthX = from.getMaxXIndex()+1;
			Size lengthY = from.getMaxYIndex()+1;
			//float stepPhysX = from.getPhysStepWidthX();
			//float stepPhysY = from.getPhysStepWidthY();
			float stepFourierX = from.getFourierStepWidthX();
			float stepFourierY = from.getFourierStepWidthY();


			
			RegularData2D newGrid(RegularData2D::IndexType(lengthX, lengthY),
														Vector2(from.getFourierSpaceMinX(), 
																		from.getFourierSpaceMinY()),
														Vector2(from.getFourierSpaceMaxX(),
																		from.getFourierSpaceMaxY()));


			// and fill it
			// AR: old version double normalization=1./(sqrt(2.*M_PI))*(stepPhysX*stepPhysY*stepPhysZ)/(pow((float)(lengthX*lengthY*lengthZ),from.getNumberOfInverseTransforms()));
			double normalization=1./(2.*M_PI)/(pow((float)(lengthX*lengthY),from.getNumberOfInverseTransforms()));
			
			Index x, y;
			signed int xp, yp;
			Vector2 r;
			typename TFFT2D<ComplexTraits>::Complex dataIn;
			typename TFFT2D<ComplexTraits>::Complex dataOut;
	
			for (Position i = 0; i < from.size(); i++)
			{
				y =  i % lengthY;
				x =  i / lengthY;
				
				xp = x;
				yp = y;

				if (xp>=lengthX/2.)
				{
					xp-=(int)lengthX;
				}
				if (yp>=lengthY/2.)
				{
					yp-=(int)lengthY;
				}

				if (x>=lengthX/2.)
				{
					x-=(int)(lengthX/2.);
				}
				else
				{
					x+=(int)(lengthX/2.);
				}

				if (y>=lengthY/2.)
				{
					y-=(int)(lengthY/2.);
				}
				else
				{
					y+=(int)(lengthY/2.);
				}


				r.set((float)xp * stepFourierX,
							(float)yp * stepFourierY);

				dataIn = from[i];
				dataOut = dataIn;

				newGrid[x + y*lengthY] = (dataOut*(typename ComplexTraits::ComplexPrecision)normalization*from.phase(r)).real();
			}

			to = newGrid;

			return to;
		}
	}
	
	
	
}
			 

#endif // BALL_MATHS_TFFT2D_H
