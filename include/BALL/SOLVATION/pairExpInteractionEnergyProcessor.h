// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// $Id: pairExpInteractionEnergyProcessor.h,v 1.19.6.2 2005/08/13 15:58:18 amoll Exp $
//

#ifndef BALL_SOLVATION_PAIREXPINTERACTIONENERGYPROCESSOR_H
#define BALL_SOLVATION_PAIREXPINTERACTIONENERGYPROCESSOR_H

#ifndef BALL_COMMON_H
# include <BALL/common.h>
#endif

#ifndef BALL_DATATYPE_OPTIONS_H
# include <BALL/DATATYPE/options.h>
#endif

#ifndef BALL_KERNEL_ATOM_H
# include <BALL/KERNEL/atom.h>
#endif

#ifndef BALL_MATHS_SURFACE_H
# include <BALL/MATHS/surface.h>
#endif

#ifndef BALL_ENERGY_ENERGYPROCESSOR_H
# include <BALL/ENERGY/energyProcessor.h>
#endif

#ifndef BALL_STRUCTURE_RDFPARAMETER_H
# include <BALL/STRUCTURE/RDFParameter.h>
#endif

#ifndef BALL_SOLVATION_PAIREXPRDFINTEGRATOR_H
# include <BALL/SOLVATION/pairExpRDFIntegrator.h>
#endif

#ifndef BALL_SOLVATION_SOLVENTDESCRIPTOR_H
# include <BALL/SOLVATION/solventDescriptor.h>
#endif

// ?????: The constants alpha, C1 and C2 i.e. the K_ij have to be
// embedded in a senseful way.

namespace BALL
{
	/** Processor for the computation of the van-derWaals interaction energy
			of a molecule.
			This processor uses a 6-exp pair potential for the calculation of
			dispersion and repulsion energies.  \par
	\ingroup Solvation		
	 */
	
	class BALL_EXPORT PairExpInteractionEnergyProcessor
		:	public EnergyProcessor
	{

		public:

		BALL_CREATE(PairExpInteractionEnergyProcessor)

		// ?????: Doku.
		enum SurfaceType
		{
			/// Unknown surface.
			SURFACE__UNKNOWN = 0,
			/// Solvent assessible surface.
			SURFACE__SAS = 1,
			/// Solvent excluding surface
			SURFACE__SES = 2,
			/// Use a surface description from a file
			SURFACE__EXTERNAL = 3
		};

		/** Symbolic names for option keys.
				This struct contains a symbolic name for each recognized key in
				PairExpInteractionEnergyProcessor::options.
		 */
		struct BALL_EXPORT Option
		{
			
			/** The verbosity level.
					Use integer values with this option.
					@see Default::VERBOSITY
					@param verbosity integer
			 */
			static const char* VERBOSITY;

			/** Constants for the pair potential.
					@see Default::ALPHA
					@see Default::C1
					@see Default::C2
					@param alpha float
					@param c1 float
					@param c2 float
			 */
			static const char* ALPHA;
			static const char* C1;
			static const char* C2;

			static const char* CLAVERIE_FILENAME;

			/** RDF option.
					This option states whether the RDF should be considered during the
					integrtion or not. Use bool values with this option.
					@see Default::USE_RDF
					@param verbosity integer
			 */
			static const char* USE_RDF;

			/** RDF file option.
					This options sets the name of the file containing the RDF
					information. Use char* values with this option.
					@see Default::RDF_FILENAME
					@param rdf_file_name char*
			 */
			static const char* RDF_FILENAME;

			/** Solvent description file option.
					This option sets the name of the file containig the solvent
					description. Use char* values with this option.
					@see Default::SOLVENT_FILENAME;
					@param solvent_file_name char*
			 */
			static const char* SOLVENT_FILENAME;

			/** The type of surface to be used.
					@see Default::SURFACE_TYPE
					@param surface_type int the type of the surface
			 */
			static const char* SURFACE_TYPE;

			/** The name of the file containing a surface definition.
					@see Default::SURFACE_FILENAME
					@param surface_filename char*
			 */
			static const char* SURFACE_FILENAME;

		};

		/** Default values for interaction energy calculations.
				These values represent the default settings for the calculations 
				of the interaction energy.
		 */
		struct BALL_EXPORT Default
		{
			/** Default verbosity level.
					@see Option::VERBOSITY
			 */
			static const Size VERBOSITY;

			/** Default pair potential constants.
					@see: Option::ALPHA
					@see: Option::C1
					@see: Option::C2
			*/
			static const double ALPHA;
			static const double C1;
			static const double C2;

			static const char* CLAVERIE_FILENAME;

			/** Default RDF setting.
					We use RDF information for the calculation of the interaction
					energy by default.
					@see Option::USE_RDF
			 */
			static const bool USE_RDF;
			static const char* RDF_FILENAME;
			static const char* SOLVENT_FILENAME;

			/** 
					@see Option::SURFACE_TYPE
			 */
			static const Size SURFACE_TYPE;

			/** 
					@see Option::SURFACE_FILENAME
			 */
			static const char* SURFACE_FILENAME;

		};

		/** @name Constructors and destructors 
		*/
		//@{

		/** Default constructor 
		*/
		PairExpInteractionEnergyProcessor() throw();

		/** Copy constructor 
		*/
		PairExpInteractionEnergyProcessor(const
				PairExpInteractionEnergyProcessor& proc) throw();

		/** Destructor 
		*/
		virtual ~PairExpInteractionEnergyProcessor() throw();

		//@}
		/** @name Assignment 
		*/
		//@{

		/** Assignment operator 
		*/
		const PairExpInteractionEnergyProcessor& operator = 
			(const PairExpInteractionEnergyProcessor& proc) throw();

		/** Clear function 
		*/
		virtual void clear() throw();
		
		//@}
		/** @name Processor functions 
		*/
		//@{

		/** 
		*/
		virtual bool finish() throw();

		//@}
		/** @name Options 
		*/
		//@{

		/** Options for the calculation of the free energy 
		*/
		Options options;

		//@}

		protected:

		/*_ Alpha constant from the Kitaygorodski Potential 
		*/
		double alpha_;

		/*_ Repulsion constant from the Kitaygorodski Potential 
		*/
		double C1_;

		/*_ Dispersion constant from the Kitaygorodski Potential 
		*/
		double C2_;

		/*_ The solvent description 
		*/
		SolventDescriptor solvent_;

		/*_ The helper class for reading rdf descriptions from an INIFile 
		*/
		RDFParameter rdf_parameter_;

		// PairExpRDFIntegrator rdf_integrator_;
	

		private:

		void computeClaverieParameters(Atom::Type solvent_type,
				Atom::Type solute_type, std::pair<float, float>& parameters)
			throw();
		void getExternalSurface_(
				std::vector< std::pair<Vector3, Surface> >& surface_map, 
				const char* surface_file) throw();

	};
   
} // namespace BALL

#endif // BALL_SOLVATION_PAIREXPINTERACTIONENERGYPROCESSOR_H
