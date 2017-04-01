// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// $Id: directory.h,v 1.27.4.2 2005/12/06 12:44:12 amoll Exp $
//

#ifndef BALL_SYSTEM_DIRECTORY_H
#define BALL_SYSTEM_DIRECTORY_H

#ifndef BALL_COMMON_H
#	include <BALL/common.h>
#endif

#ifndef BALL_CONCEPT_PROCESSOR_H
#	include <BALL/CONCEPT/processor.h>
#endif

#ifndef BALL_DATATYPE_STRING_H
#	include <BALL/DATATYPE/string.h>
#endif

#ifndef BALL_SYSTEM_FILESYSTEM_H
#	include <BALL/SYSTEM/fileSystem.h>
#endif

#ifdef BALL_HAS_DIRENT_H
#	include <dirent.h>
#endif
#ifdef BALL_HAS_UNISTD_H
#	include <unistd.h>
#endif
#ifdef BALL_HAS_SYS_STAT_H
#	include <sys/stat.h>
#endif
#ifdef BALL_HAS_DIRECT_H
#	include <direct.h>
#endif
#include <stdio.h>
#ifdef BALL_COMPILER_MSVC
#include <windows.h>
#endif

namespace BALL 
{
	/**	Directory class. 
	\ingroup System		
	*/
	class BALL_EXPORT Directory
	{
#ifndef BALL_PLATFORM_WINDOWS
# define INVALID_HANDLE_VALUE 0 
#endif

		public:

		/**	@name Constants
		*/
		//@{

		/**	The maximum length of a path.
				This constant is used for system calls that require
				a maximum length (e.g., getcwd()). Default is 8192.
		*/
		static const Size MAX_PATH_LENGTH;

		//@}
		/**	@name	Constructors and Destructors 
		*/
		//@{

		/** Default constructor.
				Constructs new Directory object.
				The directory path is set to the current working directory.
				The directory path does not have a path seperator {\em "/"} at its end.
				@return    Directory - new constructed Directory object
		*/
		Directory();

		/** Detailed constructor.
				Construct new Directory object from {\em directory_path}.
				The directory path does not have a path seperator {\em "/"} at its end.
				If the given directory does not exists, the directory path is set to an empty string.
				@param  directory_path the name of the directory to be opend
				@param  set_current true, to set the directory as the current, default = false
				@return Directory - new constructed Directory object
		*/
		Directory(const String& directory_path, bool set_current = false);

		/** Copy constructor.
				Construct new Directory object by copying <b>  directory </b>.
				@param  directory the Directory object to be copied (cloned)
				@return Directory - new constructed Directory cloned from <b>  directory </b>
		*/
		Directory(const Directory& directory);

		/** Destructor.
		*/
		~Directory();

		/** Explicit default initialization.
				Set the state to the default values.
				The path is set to an empty string.
		*/
		void clear();

		/** Explicit destructor.
		*/
		void destroy();

		//@}
		/**	@name	Assignment 
		*/
		//@{

		/** Assign the Directory with the path <tt>directory_path</tt>.
				The given directory path can be either absolute or relative. 
				If the path starts with a path seperator it is set as a absolute path.
				@param  directory_path the name of the directory to be cloned
				@param  set_current true to set the directory as the current, default = false
				@return bool true if the path could be set and is valid
		*/
		bool set(const String& directory_path, bool set_current = false);

		/** Assignment with cloning facility.
				Assign the Directory <tt> directory</tt> to <tt> *this</tt>.
				@param  directory the directory to be cloned
		*/
		void set(const Directory& directory);

		/** Assignment operator.
				Assign <b>  directory </b> to this instance.
		*/
		Directory& operator = (const Directory& directory);

		/** Copying with cloning facility.
				Copy this instance to <b>  directory </b>.
				@param directory the directory to be assigned to
		*/
		void get(Directory& directory) const;
		//@}

		/**	@name	Accessors 
		*/
		//@{

		/** Get the path of this instance.
				The directory path does not have a path seperator {\em "/"} at its end
				and is absolute. If a unvalid path was set the path is an empty string.
				@return String the name of the directory
		*/
		const String& getPath() const;

		/** Rename a given directory.
				With this method the directory associated with this object can not
				be renamed. Use renameTo instead to do so.
				@param old_path the old path
				@param new_path the new path
				@return bool  true if the directory could be renamed
		*/
		bool rename(String old_path, String new_path);

		/** Rename the directory associated with this object.
				@param new_path the new path
				@return bool  true if the directory could be renamed
		*/
		bool renameTo(String new_path);

		/** Set a directory as the current.
				@param directory_path the name of the directory
				@return bool true if the directory could be set as the current
		*/
		bool setCurrent(String directory_path);

		/** Set this directory as the current working directory.
				@return bool true if the directory could be set as the current
		*/
		bool setCurrent();

		/** Create a new directory.
				The directory is created using an absolute path, if it starts
				with a path seperator, else it is created in this directory.
				@param path the path of the new directory
				@param mode the access mode of the directory
				@return bool true if the directory could be created
		*/
		bool create(String path, const mode_t& mode = 0777);

		/** Remove a directory.
				With this method the directory associated with this object can not
				be removed. Use remove() instead to do so.
				@param old_path the path of the directory
				@return bool true if the directory could be removed
		*/
		bool remove(String old_path);

		/** Remove this directory.
				The directory this object points to is deleted and the object is cleared.
				@return bool true if the directory could be removed
		*/
		bool remove();

		/** Get the name of the first entry in the directory.
				@param entry reference to the name of the first entry
				@return bool true if an entry was found
		*/
		bool getFirstEntry(String& entry);

		/** Get the name of the next entry in the directory.
				@param entry reference to the name of the next entry
				@return bool true if an entry was found
		*/
		bool getNextEntry(String& entry);

		/** Count all items in the directory.
				@return Size the number of items (files, links, directories)
		*/
		Size countItems();

		/** Count the files in the directory.
				@return Size the number of files
		*/
		Size countFiles();

		/** Count the subdirectories in the directory.
				@return Size the number of subdirectories
		*/
		Size countDirectories();

		/** Find a file in the directory.
				The search is recursive.
				@param filename the name of the file to be searched
				@param filepath	the path of the file, if it was found
				@return bool true if the file was found
		*/
		bool find(const String& filename, String& filepath);

		//@}
		/**	@name	Predicates 
		*/
		//@{
		
		/** Test if the directory has an item.
				@param item the name of the item to look for
				@return bool true if the directory has the item
		*/
		bool has(const String& item);

		/**	Test if the directory is valid.
				The directory is valid if it exists.
				This function uses ::opendir(const char *dirname).
				@return bool true if the directory is valid
		*/
		bool isValid() const;

		/** Test if the directory is the current working directory.
				@return bool
		*/
		bool isCurrent() const;

		/** Test if the directory is empty.
				@return bool
		*/
		bool isEmpty();

		/**	Equality operator.
				@return bool, <b>true</b> if the name of both directories are equal
		*/
		bool operator == (const Directory& directory) const;

		/**	Inequality operator.
				@return bool, <b>true</b> if the name of both directories are inequal
		*/
		bool operator != (const Directory& directory) const;

		/// Get the home directory of the current user
		static String getUserHomeDir()
			throw();
		
		/// Goto the home directory of the current user
		static bool changeToUserHomeDir()
			throw();

		//@}

		private:
		
		//_switch to this dir
		void synchronize_();

		//_switch back to the working directory
		bool desynchronize_(bool result = true);
#ifdef BALL_COMPILER_MSVC
		HANDLE					dirent_;
		HANDLE					dir_;
#else
		DIR*						dir_;
		dirent*					dirent_;
#endif
		String					directory_path_;
		String					backup_path_;
	};

#	ifndef BALL_NO_INLINE_FUNCTIONS
#		include <BALL/SYSTEM/directory.iC>
#	endif
  
} // namespace BALL 

#endif // BALL_SYSTEM_DIRECTORY_H 
