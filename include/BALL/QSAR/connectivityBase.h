// -*- Mode: C++; tab-width: 2: -*-
// vi: set ts=2:
//
//

#ifndef BALL_QSAR_CONNECTIVITYBASE_H
#define BALL_QSAR_CONNECTIVITYBASE_H

#ifndef BALL_CONCEPT_MOLECULE_H
#include <BALL/KERNEL/molecule.h>
#endif

#ifndef BALL_QSAR_DESCRIPTOR_H
#include <BALL/QSAR/descriptor.h>
#endif

namespace BALL
{

	/** Generic QSAR molecular connectivity descriptors class
			\\
	*/
	class BALL_EXPORT ConnectivityBase
		:	public Descriptor
	{
		public:
				
		/** @name Constructors and Destructors
		*/
		//@{
		/** Default constructor
		*/
		ConnectivityBase();

		/** Copy constructor
		*/
		ConnectivityBase(const ConnectivityBase& cb);

		/** Named constructor
		*/
		ConnectivityBase(const String& name);

		/** Named unity constructor
		*/
		ConnectivityBase(const String& name, const String& unit);

		/** Destructor
		*/
		virtual ~ConnectivityBase();
		//@}

		/** @name Assignment
		*/
		//@{
		/** Assignment operator
		*/
		virtual ConnectivityBase& operator = (const ConnectivityBase& cb);

		protected:

		/** @name Predicates
		*/
		//@{
		bool isValid(Molecule& molecule);
		//@}
		
		/** @name Accessors
		*/
		//@{
		void calculate(Molecule& molecule);
		//@}


		private:

		/*_ @name Accessors
		*/
		//@{
		/*_ Dijkstra recursion. Performs a single source shortest path approach
				@param distances to the other nodes are stored in referenced vector<double>
				@param Atom which acts as source atom
				@param indices map, which maps for each atom a index for the dist vector
		*/
		void recursion_(std::vector<double>& dists, const Atom* source, HashMap<const Atom*, Size>& index_map);
		//@}

	};

} // namespace BALL

#endif
