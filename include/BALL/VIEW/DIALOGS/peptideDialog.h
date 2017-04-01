// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
//

#ifndef BALL_VIEW_DIALOGS_PEPTIDEDIALOG_H
#define BALL_VIEW_DIALOGS_PEPTIDEDIALOG_H

#include <BALL/VIEW/UIC/peptideDialogData.h>

#ifndef BALL_STRUCTURE_PEPTIDEBUILDER_H
#include <BALL/STRUCTURE/peptideBuilder.h>
#endif

namespace BALL
{
	namespace VIEW
	{
		/** Dialog for building peptides
				\ingroup ViewDialogs
		*/
		class BALL_VIEW_EXPORT PeptideDialog
			: public PeptideDialogData,
				public Peptides::PeptideBuilder
		{ 
			Q_OBJECT

			public:

			/// Constructor
			PeptideDialog( QWidget* parent = 0, const char* name = 0, WFlags fl = 0 );

			/// Destructor
			~PeptideDialog() {}

			///
			String getSequence()
				throw();

			///
			Protein* getProtein() 
				throw() { return protein_;}

			public slots:

			///
			void ala_pressed() {insert_('a');}

			///
			void arg_pressed() {insert_('r');}
			
			///
			void asn_pressed() {insert_('n');}
			
			///
			void asp_pressed() {insert_('d');}
			
			///
			void cys_pressed() {insert_('c');}
			
			///
			void gln_pressed() {insert_('q');}
			
			///
			void glu_pressed() {insert_('e');}
			
			///
			void gly_pressed() {insert_('g');}
			
			///
			void his_pressed() {insert_('h');}
			
			///
			void ile_pressed() {insert_('i');}
			
			///
			void leu_pressed() {insert_('l');}
			
			///
			void lys_pressed() {insert_('k');}
			
			///
			void met_pressed() {insert_('m');}
			
			///
			void phe_pressed() {insert_('f');}
			
			///
			void pro_pressed() {insert_('p');}
			
			///
			void ser_pressed() {insert_('s');}
			
			///
			void thr_pressed() {insert_('t');}
			
			///
			void trp_pressed() {insert_('w');}
			
			///
			void tyr_pressed() {insert_('y');}
			
			///
			void val_pressed() {insert_('v');}

			///
			void back_pressed();

			///
			void close_pressed();

			///
			void angle_changed();
			
			///
			void insert_seq();

			protected:
 
			void insert_(char a)
				throw();

			void truncString_(String& data)
				throw();

			bool isAminoAcid_(char a) const;

			Protein* protein_;
			static String all_amino_acids_;
		};

} } // namespaces
#endif
