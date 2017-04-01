// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// $Id: XYZFile.h,v 1.16.6.2 2005/08/12 14:03:15 amoll Exp $
//

#ifndef BALL_FORMAT_XYZFILE_H
#define BALL_FORMAT_XYZFILE_H

#ifndef BALL_SYSTEM_FILE_H
#	include <BALL/SYSTEM/file.h>
#endif

namespace BALL 
{
	class System;

	/**	XYZ file class.
			This class enables BALL to read and write XMol XYZ files.
			The XYZ format is a very simple molecular file format. It contains
			only the atom type (i.e., the element) and the cartesian coordinates
			of the structure. Bonds, atom names, and structural informations are
			not part of this file format. \par
			The first line of each XYZ file contains a single integer number: the number
			of atoms in the file. The second line is just a comment line. When reading a 
			XYZ file, BALL stores this comment as the name attribute of the system read.
			Similarly, on writing the system, it's name is written to this comment line.
			All remaining lines contain the element symbol and the three coordinates
			in free format.	 \par
			
    	\ingroup  StructureFormats
	*/
	class BALL_EXPORT XYZFile
		: public File
	{
		public:

		/**	@name	Constructors and Destructors
		*/
		//@{

		/**	Default constructor
		*/
		XYZFile();

		/** Detailed constructor.
				Create a XYZ file and open it with mode <tt>open_mode</tt> (reading is default)
				@param filename the filename
				@param open_mode the openmode - default is  \link File::IN File::IN \endlink 
		*/
		XYZFile(const String& filename, File::OpenMode open_mode = std::ios::in)
			throw(Exception::FileNotFound);

		/**	Copy constructor
		*/
		XYZFile(const XYZFile& file)
			throw(Exception::FileNotFound);

		/// Destructor
		virtual ~XYZFile()
			throw();
		
		//@}

		/**	@name Reading and Writing of Kernel Datastructures
		*/
		//@{
		
		/**	Write a system to the XYZ file
		*/
		virtual bool write(const System&	system)
			throw(File::CannotWrite);
		
		/**	Read a system from the XYZ file
		*/
		virtual void read(System&	system);

		/**	Read a system from the XYZ file
		*/
		virtual XYZFile& operator >> (System& system)
		{
			read(system);
			return *this;
		}
		
		/**	Write a system to the XYZ file
		*/
		virtual XYZFile& operator << (const System& system)
		{
			write(system);
			return *this;
		}

		const XYZFile& operator = (const XYZFile& file)
			throw()
		{
			File::operator = (file);

			return *this;
		}

		
		//@}
	};
} // namespace BALL

#endif // BALL_FORMAT_XYZFILE_H
