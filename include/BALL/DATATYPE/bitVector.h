// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// $Id: bitVector.h,v 1.46.4.3 2005/08/12 00:32:44 amoll Exp $
//

#ifndef BALL_DATATYPE_BITVECTOR_H
#define BALL_DATATYPE_BITVECTOR_H

#ifndef BALL_COMMON_H
#	include <BALL/common.h>
#endif

#ifndef BALL_CONCEPT_PERSISTENCEMANAGER_H
#	include <BALL/CONCEPT/persistenceManager.h>
#endif

#ifndef BALL_COMMON_EXCEPTION_H
#	include <BALL/COMMON/exception.h>
#endif


#include <string.h>

#define BALL_BLOCK_BITS 8
#define BALL_BLOCK_MASK (BALL_BLOCK_BITS - 1)
#define BALL_BLOCK_SHIFT              3
#define BALL_BLOCK_ALL_BITS_SET       0xFF
#define BALL_BLOCK_ALL_BITS_CLEARED   0x00


#define BALL_BLOCK_SIZE(bits) (Size)(((bits) + BALL_BLOCK_BITS - 1) >> BALL_BLOCK_SHIFT)


namespace BALL 
{
	class BitVector;

	/**	Bit Class. This class represents a bit within a BitVector.
			
			@see BitVector. 
    	\ingroup  DatatypeMiscellaneous
	*/
	class BALL_EXPORT Bit
	{
		public:

		/**	@name Exceptions
		*/
		//@{

		/** Exception thrown if a file could not be processed right.
		*/
		class BALL_EXPORT IllegalOperation
			: public Exception::GeneralException
		{
			public:
			IllegalOperation(const char* file, int line);
		};

		//@}


		/**	@name	Constructors and Destructors
		*/
		//@{

		BALL_CREATE(Bit)

		/** Default constructor 
		*/
		Bit() throw();

		/**	Copy constructor
		*/
		Bit(const Bit& bit) throw();

		/** Detailed constructor.
				For use with nonconst bitvector.
				The bitvector can be resized by accessing bits out of the bitvectors range.
				@exception NullPointer if bitvector is equal to 0
		*/
		Bit(BitVector* bitvector, Index index = 0) 
			throw(Exception::NullPointer);

		/** Detailed constructor.
				For use with const bitvector.
				@exception NullPointer if bitvector is equal to 0
				@exception IndexUnderflow if index is too small
				@exception IndexOverflow	if index is greater than the size of bitvector
		*/
		Bit(const BitVector* const bitvector, Index index = 0) 
			throw(Exception::NullPointer, Exception::IndexUnderflow, Exception::IndexOverflow);

		/** Destructor
		*/
		virtual ~Bit()	
			throw();
		//@}


		/**	@name	Converters
		*/
		//@{
		
		/** Casting operator from Bit to bool.
		*/
		operator bool() const throw(Exception::NullPointer);

		//@}
		/**	Assignment
		*/
		//@{
		
		/** Assignment operator.
				Assign the position from a Bit to this instance
		*/
		Bit& operator = (const Bit& bit) throw();

		/** Assignment operator.
				Assign a bool value to this instance.
				The bit in the bitvector is set to the given value.
				@exception IllegalOperation if instance points to a const bitvector
		*/
		Bit& operator = (const bool bit) 
			throw(Exception::NullPointer, IllegalOperation);

		/** Clear method 
		*/
		virtual void clear() throw();

		//@}
		/**	@name	 Predicates
		*/
		//@{

		/** Equality operator.
				Test if two instances have the same position in a bitvector
		*/
		bool operator == (const Bit& bit) const throw();

		/** Equality operator.
				Test if this instance has the given bool value
		*/
		bool operator == (bool bit) const throw(Exception::NullPointer);

		/** Inequality operator.
				Test if two instances point to different positions.
		*/
		bool operator != (const Bit& bit) const throw();

		/** Inequality operator.
				Test if this instance has not the given bool value
		*/
		bool operator != (bool bit) const throw(Exception::NullPointer);

		//@}

		private:

		// --- ATTRIBUTES

		BitVector* 	bitvector_;
		Index 			index_;
		bool				bitvector_muteable_;
	};


	/**	Bit vector class.
			Indices may be given as negative arguments: start from the end
			-1 therefore means the last bit.
			Some functions resize the instance if a index greater than the
			size of the instance is given.
			Other functions throw exception in this case.

    	\ingroup  DatatypeMiscellaneous
	*/
	class BALL_EXPORT BitVector
	{
		public:

		BALL_CREATE(BitVector)


		/**	@name	Type Definitions and Constants
		*/
		//@{

		///
		typedef unsigned char BlockType;
		///
		typedef std::vector<BlockType>	VectorType;
		///
		static const Size BlockSize;

		//@}
		/**	@name	Constructors and Destructors 
		*/
		//@{

		/**	Default constructor
		*/
		BitVector() throw();

		/**	Detailled constructor
		*/
		BitVector(Size size) throw(Exception::OutOfMemory);

		/** Copy constructor
		*/
		BitVector(const BitVector& bit_vector) 
			throw(Exception::OutOfMemory);

		/** Detailled constructor from an array of char.
				{\em Caveat:} the array of char has to be zero-terminated!
		*/
		BitVector(const char* bit_string)
		 throw(Exception::OutOfMemory);

		/**	Destructor
		*/
		virtual ~BitVector() throw();

		/**	Clear method.
				This method will set the size of this instance to 0.
		*/
		void clear() throw();

		//@}
		/**	@name	Assignment 
		*/
		//@{

		/// Assignment from an other BitVector instance.
		void set(const BitVector& bit_vector)
		  throw(Exception::OutOfMemory);

		/** Assignment from a char string.
				If a char is different from '0', it is interpreted as true.
				{\em Caveat:} the array of char has to be zero-terminated!
		*/
		void set(const char* bit_string)
		  throw(Exception::OutOfMemory);

		/// Assignment from an other BitVector instance.
		BitVector& operator = (const BitVector& bit_vector)
		  throw(Exception::OutOfMemory);

		/** Assignment from a char string.
				If a char is different from '0', it is interpreted as true.
		*/
		BitVector& operator = (const char *bit_string)
		  throw(Exception::OutOfMemory);

		/// Assignment to an other BitVector.
		void get(BitVector& bitvector) const
		  throw(Exception::OutOfMemory);
		//@}

		/**	@name	Accessors 
		*/
		//@{

		/** Return a sub-BitVector.
				A new BitVector is created and filled with elements of this instance.
				@param first the index of the first element to be copied
				@param last the index of the last element to be copied
				@return BitVector a partial copy of this instance
		*/
		BitVector operator () (Index first,Index last) const
		  throw(Exception::IndexUnderflow, Exception::IndexOverflow);

		/** Set the size of this instance.
				@param size the new size
				@param keep ?????
		*/
		void setSize(Size size, bool keep = true)
		  throw(Exception::OutOfMemory);

		/**	Return the size of this instance.
		*/
		Size getSize() const throw();

		/** Count a given bool value in this instance.
				@param value the value to be counted
				@return Size the number of the given value
		*/
		Size countValue(bool value) const throw();

		/** Return a muteable pointer to the values of this instance.
				@return unsigned char*	a muteable pointer to the values of this instance
		*/
		VectorType& getBitSet() throw();

		/** Return a constant pointer to the values of this instance.
				@return unsigned char*	a constant pointer to the values of this instance
		*/
		const VectorType& getBitSet() const throw();

		/**	Random access to the components.
				If the given index is greater than the size of this instance, this BitVector
				is increased to the given index.
				@return Bit a Bit pointing to the given element of this instance
		*/
		Bit operator []	(Index index)
		  throw(Exception::OutOfMemory);

		/**	Constant random access to the components.
				If the given index is greater than the size of this instance an exception is thrown.
				@return bool the value of the given index
		*/
		bool operator [] (Index index) const
			throw(Exception::IndexUnderflow, Exception::IndexOverflow);

		/** Set the given Bit.
				The element at the given position is set to the given value.
				If the given index is greater than the size of this instance, this BitVector
				is increased to the given index.
				@param index the index of the element
				@param value the value to be set
		*/
		void setBit(Index index, bool value = true)
		  throw(Exception::IndexUnderflow, Exception::OutOfMemory);

		/** Get the value of an element.
				If the given index is greater than the size of this instance, this BitVector
				is increased to the given index.
				@param index the index of the element
				@return bool the value of the element
		*/
		bool getBit(Index index)
		  throw(Exception::IndexUnderflow, Exception::OutOfMemory);

		/** Get the value of an element (const version).
				If the given index is greater than the size of this instance an exception is thrown.
				No resizing is done.
				@param index the index of the element
				@return bool the value of the element
		*/
		bool getBit(Index index) const
			throw(Exception::IndexUnderflow, Exception::IndexOverflow);

		/** Toggle the value of an element.
				If the given index is greater than the size of this instance an exception is thrown.
				@param index the index of the element
		*/
		void toggleBit(Index index)
		  throw(Exception::IndexUnderflow, Exception::OutOfMemory);

		/** Fill a part of this instance with a given value.
				If a given index is greater than the size of this instance an exception is thrown.
				@param value the value to be set
				@param first the starting position
				@param last the last element to be filled
		*/
		void fill(bool value = true, Index first = 0 , Index last = -1)
		  throw(Exception::IndexUnderflow, Exception::OutOfMemory);

		/** Toggle the values of a part of this instance.
				If a given index is greater than the size of this instance an exception is thrown.
				@param value the value to be set
				@param first the starting position
				@param last the last element to be filled
		*/
		void toggle(Index first = 0, Index last = -1)
			throw(Exception::IndexUnderflow, Exception::OutOfMemory);

		/** Set a unsigned char as the bit pattern.
				For example: 22 => 00010110.
				@param bit_pattern the new pattern.
		*/
		void setUnsignedChar(unsigned char bit_pattern) throw();

		/** Get the bit pattern as unsigned char.
				For example: 00010110 => 22.
				@return unsigned char the bit-pattern
		*/
		unsigned char getUnsignedChar() const throw();

		/** Set a unsigned short as the bit pattern.
				@see setUnsignedChar
		*/
		void setUnsignedShort(unsigned short bit_pattern) throw();

		/** Get the bit pattern as unsigned short.
				@see getUnsignedChar
		*/
		unsigned short getUnsignedShort() const throw();

		/** Set a unsigned int as the bit pattern.
				@see setUnsignedChar
		*/
		void setUnsignedInt(unsigned int bit_pattern) throw();

		/** Get the bit pattern as unsigned int.
				@see getUnsignedChar
		*/
		unsigned int getUnsignedInt() const throw();

		/** Set a unsigned long as the bit pattern.
				@see setUnsignedChar
		*/
		void setUnsignedLong(unsigned long bit_pattern) throw();

		/** Get the bit pattern as unsigned long.
				@see getUnsignedChar
		*/
		unsigned long getUnsignedLong() const throw();

		/** Compute {\em this or bit_vector}.
				The result is saved in this instance.
		*/
		void bitwiseOr(const BitVector& bit_vector) 
			throw(Exception::OutOfMemory);

		/** Compute {\em this xor bit_vector}.
				The result is saved in this instance.
		*/
		void bitwiseXor(const BitVector& bit_vector) 
			throw(Exception::OutOfMemory);

		/** Compute {\em this and bit_vector}.
				The result is saved in this instance.
		*/
		void bitwiseAnd(const BitVector& bit_vector) 
			throw(Exception::OutOfMemory);

		/** Or Operator.
				Creates a new BitVector object and fills it with the result of
				{\em this or bit_Vector}.
		*/
		BitVector operator | (const BitVector& bit_vector) 
			throw(Exception::OutOfMemory);

		/** Compute {\em this or bit_vector}.
				The result is saved in this instance.
		*/
		BitVector& operator |= (const BitVector& bit_vector) 
			throw(Exception::OutOfMemory);

		/** And Operator.
				Creates a new BitVector object and fills it with the result of
				{\em this and bit_Vector}.
		*/
		BitVector operator & (const BitVector& bit_vector) 
			throw(Exception::OutOfMemory);

		/** Compute {\em this and bit_vector}.
				The result is saved in this instance.
		*/
		BitVector& operator &= (const BitVector& bit_vector) 
			throw(Exception::OutOfMemory);

		/** Xor Operator.
				Creates a new BitVector object and fills it with the result of
				{\em this Xor bit_vector}.
		*/
		BitVector operator ^ (const BitVector& bit_vector) 
			throw(Exception::OutOfMemory);

		/** Compute {\em this xor bit_vector}.
				The result is saved in this instance.
		*/
		BitVector& operator ^= (const BitVector& bit_vector) 
			throw(Exception::OutOfMemory);

		/** Negate Operator.
				Creates a new BitVector object and fills it with the negate 
				result of this instance.
		*/
		BitVector operator ~ () 
			throw(Exception::OutOfMemory);

		//@}
		/**	@name	Predicates
		*/
		//@{

		/// Equality operator
		bool operator == (const BitVector& bit_vector) const throw();

		/// Inequality operator
		bool operator != (const BitVector& bit_vector) const throw();

		/** Test if any bit in a given range has the given value.
				@param value the value to look for
				@param first the index to start searching
				@param last the index to stop searching
		*/
		bool isAnyBit(bool value, Index first = 0, Index last = -1) const
		  throw(Exception::IndexUnderflow, Exception::IndexOverflow);

		/** Test if every bit in a given range has the given value
				@param value the value to look for
				@param first the index to start searching
				@param last the index to stop searching
		*/
		bool isEveryBit(bool value, Index first = 0, Index last = -1) const
		  throw(Exception::IndexUnderflow, Exception::IndexOverflow);

		//@}
		/**	@name	Debugging and Diagnostics 
		*/
		//@{

		/** Test if this instance is valid.
		*/
		bool isValid() const throw();

		//@}
		/**	@name	Storers 
		*/
		//@{

		/**	Input operator.
				Reads the values of type <b>  bool </b> from an istream.
		*/
		BALL_EXPORT friend std::istream& operator >> (std::istream& s, BitVector& bit_vector)
			throw(Exception::OutOfMemory);

		/**	Output operator.
				Writes the values of type <b>  bool </b> to an ostream.
		*/
		BALL_EXPORT friend std::ostream& operator << (std::ostream& s, const BitVector& bit_vector)
			throw();

		/**	Read the values of of type <b>  bool </b> from an istream.
		*/
		virtual void read(std::istream& s)
			throw(Exception::OutOfMemory);

		/**	Write the values of of type <b>  bool </b> to an ostream.
		*/
		virtual void write(std::ostream& s) const throw();

		/**	Storable write method.
		*/
		virtual void write(PersistenceManager& pm) const throw();

		/**	Storable read method.
		*/
		virtual bool read(PersistenceManager& pm)
		  throw(Exception::OutOfMemory);

		//@}

		protected:
	
		void validateIndex_(Index& index)
			throw(Exception::IndexUnderflow, Exception::OutOfMemory);

		void validateIndex_(Index& index) const
		  throw(Exception::IndexUnderflow, Exception::IndexOverflow);

		void validateRange_(Index& first, Index& last) const
		  throw(Exception::IndexUnderflow, Exception::IndexOverflow);


		private:
	
		Index block_(Index index)
			throw(Exception::IndexUnderflow, Exception::OutOfMemory);

		Index block_(Index index) const
			throw(Exception::IndexUnderflow, Exception::IndexOverflow);

		BlockType mask_(Index index) const
			throw();

		// --- ATTRIBUTES

		Size		size_;
		VectorType	bitset_;
	};

#	ifndef BALL_NO_INLINE_FUNCTIONS
#		include <BALL/DATATYPE/bitVector.iC>
#	endif
} //namespace BALL


#endif // BALL_DATATYPE_BITVECTOR_H
