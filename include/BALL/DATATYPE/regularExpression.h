// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// $Id: regularExpression.h,v 1.2.6.4 2005/11/04 08:07:04 oliver Exp $
//

#ifndef BALL_DATATYPE_REGULAREXPRESSION_H
#define BALL_DATATYPE_REGULAREXPRESSION_H

#ifndef BALL_COMMON_H
#	include <BALL/common.h>
#endif

// for SGI CC: turn of warnings #1174: function "regcomp" was declared but never referenced
#if !defined(__GNUC__) && !defined(__KAI__) && defined(IRIX)
#	pragma set woff 1174
#endif

#ifdef BALL_HAS_SYS_TYPES_H
#	include <sys/types.h>
#endif

#ifdef BALL_HAS_REGEX_H
	// make sure __STDC__ is set when including the
  // GNU regex headers under Windows
# ifdef BALL_COMPILER_MSVC
#		define __STDC__ 1
		extern "C" {
#	endif
#			include <regex.h>
#	ifdef BALL_COMPILER_MSVC
		}
#	endif
#endif

#if !defined(__GNUC__) && !defined(__KAI__) && defined(IRIX)
#	pragma reset woff 1174
#endif

#ifndef BALL_DATATYPE_STRING_H
#	include <BALL/DATATYPE/string.h>
#endif

#define BALL_REGULAR_EXPRESSION_DEFAULT_PATTERN   ""

namespace BALL 
{
	/**	Regular expression class
    	\ingroup  DatatypeMiscellaneous
	*/
	class BALL_EXPORT RegularExpression
	{
		public:

		BALL_CREATE(RegularExpression)

		/**	@name	String constants 
		*/
		//@{
			
		/**	alphabetic letters.
		*/
		static const String ALPHA;

		/**	alphanumeric letters
		*/
		static const String ALPHANUMERIC;

		/**	real numbers
		*/
		static const String REAL;

		/** identifier characters
		*/
		static const String IDENTIFIER;

		/** integer characters
		*/
		static const String INTEGER;

		/** hexadecimal and integer characters
		*/
		static const String HEXADECIMAL_INTEGER;

		/** lowercase letters
		*/
		static const String LOWERCASE;

		/** non alphabetic characters
		*/
		static const String NON_ALPHA;

		/** non alphanumeric characters
		*/
		static const String NON_ALPHANUMERIC;

		/** non numeric characters
		*/
		static const String NON_NUMERIC;

		/** non whitespace charcters
		*/
		static const String NON_WHITESPACE;

		/**	uppercase letters.
		*/
		static const String UPPERCASE;

		/** whitespace characters.
		*/
		static const String WHITESPACE;
		//@}

		/**	@name	Constructors and Destructors 
		*/
		//@{
			
		///	Default constructor
		RegularExpression()	throw();

		///	Copy constructor
		RegularExpression(const RegularExpression& regular_expression) throw();

		///
		RegularExpression(const String& pattern, bool wildcard_pattern = false) throw();

		///	Destructor
		virtual ~RegularExpression() throw();

		///	Reset the object attributes to their default values
		virtual void clear() throw();

		///	Destroy the instance
		void destroy() throw();
		//@}

		/**	@name	Assignment
		*/
		//@{

		/// Assignment operator
		RegularExpression& operator = (const RegularExpression& expression) throw();

		///	Assign from another instance
		void set(const RegularExpression& regular_expression) throw();

		///	Assign from a string
		void set(const String& pattern, bool wildcard_pattern = false) throw();

		///	Assign to another instance
		void get(RegularExpression& regular_expression) const throw();
		//@}

		/**	@name	Accessors
		*/
		//@{

		///	Get the expression pattern.
		const String& getPattern() const throw();

		///	Count subexpressions.
		Size countSubexpressions() const throw();

		/**	Match a text with a given pattern.
				@param text to process
				@param pattern to compare with
				@param compile_flags ?????
				@param execute_flags ?????
				@exception NullPointer if <tt>text</tt> or <tt>pattern</tt> are NULL
		*/
		static bool match(const char* text, const char* pattern,
											int compile_flags = 0 | REG_EXTENDED | REG_NOSUB, int execute_flags = 0)
			throw(Exception::NullPointer);
		
		/**	Match a text with this regular expression.
				@param text to process
				@param from index in the string to start the matching
				@param execute_flags ?????
		*/
		bool match(const String& text, Index from = 0, int execute_flags = 0) const
			throw(Exception::NullPointer, Exception::IndexUnderflow, Exception::IndexOverflow);
	
		/**	Match a substring with this regular expression.
				@param text to process
				@param from index in the substring to start the matching
				@param execute_flags ?????
		*/
		bool match(const Substring& text, Index from = 0, int execute_flags = 0) const
			throw(Substring::InvalidSubstring, Exception::IndexUnderflow, Exception::IndexOverflow);
	
		/**	Match a C-String with this regular expression.
				@param text to process
				@param execute_flags ?????
		*/
		bool match(const char* text, int execute_flags = 0) const
			throw(Exception::NullPointer);
	
		/**	Find this expression in a string
				@param from index in the string to start the matching
				@param found the result is stored as a substring here
				@param execute_flags ?????
		*/
		bool find(const String& text, Substring& found,
							Index from = 0, int execute_flags = 0) const
							throw(Exception::IndexUnderflow, Exception::IndexOverflow);
			
		/**	Find this expression in a string
				@param text to process
				@param subexpressions the results are stored here
				@param from index in the string to start the matching		
				@param execute_flags ?????
		*/
		bool find(const String& text, vector<Substring>& subexpressions,
							Index from = 0, int execute_flags = 0) const
							throw(Exception::IndexUnderflow, Exception::IndexOverflow);
						
		//@}
		/**	@name	Predicates
		*/
		//@{

		/**	Test if expression is empty.
		*/
		bool isEmpty() const throw();

		/**	Equality operator
		*/
		bool operator == (const RegularExpression& regular_expression) const throw();

		/**	Inequality operator
		*/
		bool operator != (const RegularExpression& regular_expression) const throw();

		/**	Less operator
		*/
		bool operator < (const RegularExpression& regular_expression) const throw();

		/**	Less or equal operator
		*/
		bool operator <= (const RegularExpression& regular_expression) const throw();

		/**	Greater or equal operator
		*/
		bool operator >= (const RegularExpression& regular_expression) const throw();

		/**	Greater operator
		*/
		bool operator > (const RegularExpression& regular_expression) const throw();

		//@}
		/**	@name	Debugging and Diagnostics
		*/
		//@{

		/**	Test if instance is valid.
		*/
		virtual bool isValid() const throw();

		/**	Dump this instance to an ostream
				depth is normaly just used for internal use.
				@param s the ostream, default is the standard output
				@param depth the indentation depth of the output
		*/
		virtual void dump(std::ostream& s = std::cout, Size depth = 0) const
			throw();

		//@}
		/**	@name Storers
		*/
		//@{

		/**	output operator
		*/
		BALL_EXPORT
		friend std::ostream& operator << (std::ostream& s, const RegularExpression& regular_expression)
			throw();

		/**	input operator
		*/
		BALL_EXPORT
		friend std::istream& operator >> (std::istream& s, RegularExpression& regular_expression)
			throw();

		//@}
		
		private:

		void compilePattern_()
			throw();

		void toExtendedRegularExpression_() 
			throw();

			regex_t regex_;
		
		String 	pattern_;
		bool 		valid_pattern_;
	};

#	ifndef BALL_NO_INLINE_FUNCTIONS
#		include <BALL/DATATYPE/regularExpression.iC>
#	endif
} // namespace BALL

#endif // BALL_DATATYPE_REGULAREXPRESSION_H
