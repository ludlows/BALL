// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// $Id: file.h,v 1.65.4.7 2005/11/12 14:16:29 oliver Exp $
//

#ifndef BALL_SYSTEM_FILE_H
#define BALL_SYSTEM_FILE_H

#ifndef BALL_DATATYPE_REGULAREXPRESSION_H
#	include <BALL/DATATYPE/regularExpression.h>
#endif

#ifndef BALL_DATATYPE_STRING_H
#	include <BALL/DATATYPE/string.h>
#endif

#ifndef BALL_SYSTEM_FILESYSTEM_H
#	include <BALL/SYSTEM/fileSystem.h>
#endif

#include <stdlib.h>                   // 'getenv'
#include <sys/stat.h>         // 'stat', 'lstat'
#include <stdio.h>                    // 'rename'

#ifdef BALL_COMPILER_MSVC
#	define S_ISREG _S_ISREG
#	define S_ISDIR _S_ISDIR
#	define S_ISCHR _S_ISCHR
#	define S_ISBLK _S_ISBLK
#	define S_ISFIFO _S_ISFIFO
#	define access _access
#endif
#include <algorithm>
#include <iostream>
#include <fstream>
#include <sys/types.h>
#include <map>

#ifdef BALL_HAS_UNISTD_H
#	include <unistd.h>			// 'access', 'rename', 'truncate'
#endif

#ifdef BALL_COMPILER_MSVC
#	include <fcntl.h>
#	include <io.h>
	// Define the missing symbols from <unistd.h>,
	// which M$, in its infinite wisdom, was unable to provide.
#	define F_OK 0
#	define W_OK 2
#	define R_OK 4
#	ifdef IN
#		undef IN
#	endif
#	ifdef	OUT
#		undef OUT
#	endif
#endif


namespace BALL 
{
	/**	This class handles automatic file transformation methods.
		   \link File File \endlink  provides the ability to transform files on the fly using 
			predefined transformation commands  (e.g. unix-style filters). For example, compressed 
      files can be automatically decompressed by calling the unic <tt>compress</tt> command.
			The respective commands are selectedvia a suitable regular expression, usually
      matching the file suffix. A frequent application for this transformation is the
			compressed storage of PDB files in the unix compressed format (<tt>*.Z</tt>).
			Transformation manager basically contains a map consisting of two strings.
			Using  \link findTransformationCommand findTransformationCommand \endlink ,  
			\link File File \endlink  can determine whether there
			is a suitable transformation command available for a given file name.
			User-defined transformation may be defined at any time using the 
			\link TransformationManager::registerTransformation registerTransformation \endlink  method of the static 
			instance of  \link TransformationManager TransformationManager \endlink 
			accessible through  \link File::getTransformationManager File::getTransformationManager \endlink .
			\ingroup System
	*/
	class BALL_EXPORT TransformationManager
	{	
		public:
			
		/**	Constructors and Destructors
		*/
		//@{
		
		/// Default constructor
		TransformationManager();
		
		/// Destructor
		~TransformationManager();
		
		//@}
		/**	Accessors
		*/
		//@{
		
		/// Add a new transformation 
		void registerTransformation(const String& pattern, const String& command);
		
		/// Delete a transformation
		void unregisterTransformation(const String& pattern);
		
		/// Find a transformation matching a given file name
		String findTransformation(const String& name) const;

		/** Apply a suitable transformation to the string.
				This first calls  \link findTransformation findTransformation \endlink  to determine the
				command string that should be applied. The <tt>name</tt> argument
				is then replaced with the contents of the matching command in the 
				TransformationManager's map. The following rules apply (in that order):

					- <tt>%s</tt> is replaced by the full content of <tt>name</tt>
					- <tt>%f</tt> is replaced by the full content of <tt>name</tt>, 
						without any file type suffix (i.e. anything after the last dot in the filename is removed)
					- <tt>%f[suffix]</tt> is replaced by the previous content of <tt>name</tt> without the <tt>suffix</tt>
					- <tt>%b</tt> and <tt>%b[suffix]</tt> like <tt>%f</tt> and <tt> %f[suffix]</tt>, except that the
								path is removed as well, so it is only the <b>  base name </b> of the file
					- <tt>%p</tt> the path to the file
					- <tt>%t</tt> a temporary file name (all occurences of <tt>%t</tt> are replace 
							with the same file name for the same invocation of  \link transform transform \endlink,
							but different file names on subsequent invocations)
		*/
		String transform(const String& name);
		//@}

		protected:
		
		/// The map containing all transformation methods
		std::map<String, String>	transformation_methods_;
	};
		
	/**	File Class.	
			\ingroup System		
	*/
	class BALL_EXPORT File
		: public std::fstream
	{
		public:

		/**	Exception CannotWrite
				A given file could not be written, either because its not open or it has a wrong open mode.
		*/
		class BALL_EXPORT CannotWrite
			: public Exception::GeneralException
		{
			public:
			CannotWrite(const char* file, int line, const String& filename)
				throw();

			~CannotWrite()
				throw();
			String getFilename() const
				throw();

			protected:
			std::string filename_;
		};

		/**	@name	Type definitions
		*/
		//@{
			
		/**	File open modes.
				This type is used to describe the standard openmodes for files
				as described in Section 27.4.2.1.4 of the ANSI C++ standard.
		*/
		typedef std::ios::openmode OpenMode;			

		//@}

		/**	@name	Constants
		*/
		//@{
		/// Open for input (default)
		static const OpenMode MODE_IN;

		/// Open for output
		static const OpenMode MODE_OUT;

		/// Append. Seek to end before each write operation
		static const OpenMode MODE_APP;

		/// Binary mode
		static const OpenMode MODE_BINARY;

		/// Seek to end directly after opening.
		static const OpenMode MODE_ATE;

		/// Truncate an existing file.
		static const OpenMode MODE_TRUNC;
		//@}

		/**	@name	Enums
		*/
		//@{

		/**	Transformation types for file.
				This enum defines some possible types for on-the-fly file transformation.
		*/
		enum Transformation
		{
			///
			TRANSFORMATION__EXEC = 1,
			///
			TRANSFORMATION__FILTER = 2,
			///
			TRANSFORMATION__URL = 3
		};

		/** Filetype
		*/
		enum Type
		{
			///
			TYPE__UNKNOWN            = 0,
			///
			TYPE__DIRECTORY          = 1,
			///
			TYPE__CHAR_SPECIAL_FILE  = 2,
			///
			TYPE__BLOCK_SPECIAL_FILE = 3,
			///
			TYPE__REGULAR_FILE       = 4,
			///
			TYPE__SYMBOLIC_LINK      = 5,
			///
			TYPE__SOCKET             = 6,
			///
			TYPE__FIFO_SPECIAL_FILE  = 7
		};

		/// Prefix for filenames that are created through the execution of commands "exec:"
		static const String TRANSFORMATION_EXEC_PREFIX;

		/// Prefix for files (to mimick URL-like behavior) "file:"
		static const String TRANSFORMATION_FILE_PREFIX;

		/// Prefix for FTP-transfers "ftp://"
		static const String TRANSFORMATION_FTP_PREFIX;

		/// Prefix for HTTP-transfer "http://"
		static const String TRANSFORMATION_HTTP_PREFIX;

		//@}
		/**	@name	Constructors and Destructors
		*/
		//@{

		/** Default constructor.
		*/
		File()
			throw();

		/** Construct new File object from the file <b>  name </b> and open the file.
				@param  name the name of the file to be opend
				@param  open_mode the openmode to be used
				@see    open
		*/
		File(const String& name, OpenMode open_mode = std::ios::in)
			throw(Exception::FileNotFound);

		/** Copy constructor.
				The file is not opend.
				@param  file the File object to be copied (cloned)
				@see    open
		*/
		File(const File& file)
			throw(Exception::FileNotFound);

		/** Destructor.
				The file is closed.
		*/
		virtual ~File()
			throw();

		/** Clear the File object.
		*/
		virtual void clear() throw();
		//@}

		/**	@name	Assignment 
		*/
		//@{

		/** Assignment operator.
				Assign the filename from <b>  file </b>.
				The file is not opend.
		*/
		const File& operator = (const File& file)
			throw();

		//@}

		/**	@name	Accessors 
		*/
		//@{

		/**	Open a given file.
				The standard constructor uses this method.
				@param name the name of the file
				@param open_mode the open mode, default is IN
				@return bool true if the file could be opened
		*/
		bool open(const String& name, File::OpenMode open_mode = std::ios::in)
			throw (Exception::FileNotFound);

		/**	Reopen the file.
				The file is closed and reopend.
				@return bool true if the file could be reopend
		*/
		bool reopen()
			throw (Exception::FileNotFound);

		/**	Reopen the file with a different mode.
				The file is closed and reopend.
				@param open_mode the new mode
				@return bool true if the file could be reopend
		*/
		bool reopen(File::OpenMode open_mode)
			throw (Exception::FileNotFound);

		/**	Close the file.
		*/
		void close()
			throw();

		/**	Return the name of the file.
				@return String the name of the file
		*/
		const String& getName()
			const	throw();
		
		/** Close the file and point to an other file.
				@param  name the new file
		*/
		void setName(const String& name)
			throw();

		/**
		*/
		const String& getOriginalName() const
			throw();

		/**	Return the size of the file.
				If the file does not exist, 0 is returned.
				@return Size the size of the file
		*/
		Size getSize()
			throw(Exception::FileNotFound);

		/**	Return the size of a given file.
				@return Size the size of the file
		*/
		static Size getSize(String name)
			throw (Exception::FileNotFound);

		/** Return the open mode.
				Default is IN.
				@return int the open mode
		*/
		File::OpenMode getOpenMode()
			const	throw();
		
		/**	Return the filetype of a given file.
				@param name the name of the file.
				@param trace_link true to follow links
				@return Type the filetype
		*/
		static Type getType(String name, bool trace_link)
			throw(Exception::FileNotFound);

		/**	Return the filetype.
				@param trace_link true to follow links
				@return Type the filetype
		*/
		Type getType(bool trace_link)
			const	throw(Exception::FileNotFound);

		/**	Copy a given file to a given destination.
				If a file with the destination name exists already, nothing happens.
				@param source_name the name of the source file
				@param destination_name the name of the destination file
				@param buffer_size the buffer size to use while copying
				@return true if copying was successfull
		*/
		static bool copy(String source_name, String destination_name, Size buffer_size = 4096)
			throw(Exception::FileNotFound);

		/**	Copy the file to a given destination.
				If a file with the destination name exists already, nothing happens.
				@param destination_name the name of the destination file
				@param buffer_size the buffer size to use while copying
				@return true if copying was successfull
		*/
		bool copyTo(const String& destination_name, Size buffer_size = 4096)
			throw(Exception::FileNotFound);

		/**	Move a given file to a given destination.
				If a file with the destination name exists already, nothing happens.
				@param source_name the name of the source file
				@param destination_name the name of the destination file
				@return true if copying was successfull
		*/
		static bool move(const String& source_name, const String& destination_name)
			throw(Exception::FileNotFound);

		/**	Move the file to a given destination.
				If a file with the destination name exists already, nothing happens.
				@param destination_name the name of the destination file
				@return true if copying was successfull
		*/
		bool moveTo(const String& destination_name)
			throw(Exception::FileNotFound);

		/**	Remove the given file.
				@param name the name of the file to be removed
				@return bool true if the file could be removed
		*/
		static bool remove(String name)
			throw();

		/**	Remove the file.
				@return bool true if the file could be removed
		*/
		bool remove()
			throw();

		/**	Rename a file.
				@param old_path the path and name of the file to be renamed
				@param new_path the new path and name of the file
				@return bool true if the file could be renamed
		*/
		static bool rename(String old_path, String new_path)
			throw (Exception::FileNotFound);

		/**	Rename the file to a given name.
				If a file with the destination name exists already, nothing happens.
				@param new_path the new path and name of the file
				@return bool true if the file could be renamed
		*/
		bool renameTo(const String& new_path)
			throw (Exception::FileNotFound);

		/**	Truncate a given file.
				@param path the path to the file
				@param size the new size of the file
				@return bool true if the file could be truncated
		*/
		static bool truncate(String path, Size size = 0)
			throw (Exception::FileNotFound);

		/**	Truncate the file.
				@param size the new size of the file
				@return bool true if the file could be truncated
		*/
		bool truncate(Size size = 0)
			throw (Exception::FileNotFound);
			
		/**	Create a temporary filename.
				This method creates strings, starting at _AAAAAAA.TMP and tries if a 
				file with this name exists. If not, the string is returned. If a file
				with this name exists, it continues to create names up to _ZZZZZZZ.TMP.
				@param temporary reference to the temporary filename
				@return bool true if a temporary filename could be found
		*/
		static bool createTemporaryFilename(String& temporary)
			throw();

    /** Return the stream associated with this file.
				Implemented just for convenience.
        @return std::fstream the stream
    */
    std::fstream& getFileStream();

		//@}
		/**@name On-the-fly file transformation
				@see TransformationManager
		*/
		//@{

		/**	Mutable access the TransformationManager.
				\link File File \endlink  defines a static instance of  
				\link TransformationManager TransformationManager \endlink  to
				handle on-the-fly conversions of files (e.g. compression, charset conversion, etc.).
		*/
		TransformationManager& getTransformationManager();

		/**	Constant access to the TransformationManager.
				\link File File \endlink  defines a static instance of  
				\link TransformationManager TransformationManager \endlink  to
				handle on-the-fly conversions of files (e.g. compression, charset conversion, etc.).
		*/
		const TransformationManager& getTransformationManager() const;
		
		/**
		*/
		static void enableTransformation(Transformation transformation);

		/**
		*/
		static void disableTransformation(Transformation transformation);

		/**
		*/
		static bool isTransformationEnabled(Transformation transformation);

		/**	
		*/
		static void registerTransformation(const String& pattern, const String& exec);

		/**	
		*/
		static void unregisterTransformation(const String& pattern);

		//@}
		/**	@name Predicates 
		*/
		//@{

		/**	Equality comparison operator.
				Two File objects are equal if they point to the same canonized filename.
		*/
		bool operator == (const File& file)
			const	throw();
		
		/**	Inequality comparison operator.
				Two File objects are inequal if they point not to the same canonized filename.
		*/
		bool operator != (const File& file)
			const	throw();

		/**	Test if the file is opend.
				The standard constructor opens the file.
				@return bool true if the file is closed
		*/
		bool isOpen()
			const	throw();

		/**	Test if the file is closed.
				The standard constructor opens the file.
				@return bool true if the file is closed
		*/
		bool isClosed()
			const	throw();

		/**	Test if a given file can be accessed.
				@param name the name of the file to be tested
				@return bool true if the file can be accessed
		*/
		static bool isAccessible(String name)
			throw ();

		/**	Test if the file can be accessed.
				@return bool true if the file can be accessed
		*/
		bool isAccessible()
			const throw (Exception::FileNotFound);

		/**	Test if the path of the file is canonized.
				The path is	compared before and after call of 
				FileSystem::canonizePath(canonized_name).
				@see FileSystem::canonizePath
				@return bool true if the path is canonized.
		*/
		bool isCanonized()
			const throw (Exception::FileNotFound);
	
		/**	Test if a given file is readable.
				@param name the name of the file
				@return true if the file can be read
		*/
		static bool isReadable(String name)
			throw (Exception::FileNotFound);

		/**	Test if the file is readable.
				@return true if the file can be read
		*/
		bool isReadable()
			const throw (Exception::FileNotFound);

		/**	Test if a given file is writeable.
				@param name the name of the file
				@return true if the file is writeable
		*/
		static bool isWritable(String name)
			throw (Exception::FileNotFound);

		/**	Test if the file is writeable.
				@return true if the file is writeable
		*/
		bool isWritable()
			const throw (Exception::FileNotFound);

		/**	Test if a given file is executable.
				@param name the name of the file
				@return true if the file is executable
		*/
		static bool isExecutable(String name)
			throw (Exception::FileNotFound);

		/**	Test if the file is executable.
				@return true if the file is executable
		*/
		bool isExecutable()
			const throw (Exception::FileNotFound);
 
		//@}
		/**	@name	Debugging and Diagnostics 
		*/
		//@{

		/**	Test if the file is valid.
				If the filename was not set, false is returned.
				This function uses std::fstream::good().
				@return bool true if the file is valid
		*/
		bool isValid()
			const	throw();

		//@}

		protected:

		String		name_;
		String		original_name_;
		OpenMode	open_mode_;
		bool			is_open_;
		bool			is_temporary_;

		static TransformationManager	transformation_manager_;
		static Size										transformation_methods_;
	};


	/** Helper class for data conversion.	
			BinaryFileAdaptors are used to read and write binary data from and to
			streams. This is done by reading the member <tt>data</tt> as a byte stream
			through an explicit cast and utilizing the stream read() and write() 
			functions. \par
			<b>Caveat:</b> This concept relies on the C++ memory layout and thus 
			is highly non-portable!
			\par
			\ingroup System		
	*/
	template <typename T>
	class BinaryFileAdaptor
	{

		public:

		/// @name Constructors and destructor
		//@{

		/// Default constructor
		BinaryFileAdaptor()
			throw();

		/// Detailed constructor
		BinaryFileAdaptor(const T& data)
			throw();

		//@}
		///@name Accessors
		//@{

		/** Sets the member <tt>data</tt> to the desired value.
				@param data data of type T
		*/
		void setData(const T& data)
			throw();

		/** Returns a const reference to the data stored in the adaptor
		*/
		const T& getData() const
			throw();

		/** Returns a mutable reference to the data stored in the adaptor
				*/
		T& getData()
			throw();

		//@}

		protected:

		//_ The member data.
		T data_;
	};

	template <typename T>
	BALL_INLINE
	BinaryFileAdaptor<T>::BinaryFileAdaptor()
		throw()
		: data_()
	{
	}

	template <typename T>
	BALL_INLINE
	BinaryFileAdaptor<T>::BinaryFileAdaptor(const T& data)
		throw()
		: data_(data)
	{
	}

	template <typename T>
	BALL_INLINE
	void BinaryFileAdaptor<T>::setData(const T& data)
		throw()
	{
		data_ = data;
	}

	template <typename T>
	BALL_INLINE
	const T& BinaryFileAdaptor<T>::getData() const 
		throw()
	{
		return data_;
	}

	template <typename T>
	BALL_INLINE
	T& BinaryFileAdaptor<T>::getData()
		throw()
	{
		return data_;
	}

	/// Output stream for BinaryFileAdaptors
	template <typename T>
	std::ostream& operator << (std::ostream& os, const BinaryFileAdaptor<T>& data)
	{
		os.write(reinterpret_cast<const char*>(&data.getData()), sizeof(T));
		return os;
	}

	/// Input stream for BinaryFileAdaptors
	template <typename T>
	std::istream& operator >> (std::istream& is, BinaryFileAdaptor<T>& data)
	{
		is.read(reinterpret_cast<char*>(&data.getData()), sizeof(T));
		return is;
	}


	/** Coping with endianness. This function swaps the bytes of a variable
			of type T if this type is of size 2n.
	*/
	template <typename T>
	BALL_INLINE
	void swapBytes(T& t)
	{
		if (sizeof(T) % 2 != 0)
		{
			Log.error() << "Cannot swap types of uneven size." << std::endl;
			return;
		}

		char* tmp = reinterpret_cast<char*>(&t);
		std::reverse(tmp, tmp + sizeof(T));
	}
	

#	ifndef BALL_NO_INLINE_FUNCTIONS
#		include <BALL/SYSTEM/file.iC>
#	endif
  
} // namespace BALL

#endif // BALL_SYSTEM_FILE_H
