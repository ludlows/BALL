// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// $Id: label.h,v 1.10.6.4 2005/09/01 22:17:59 amoll Exp $
//

#ifndef BALL_VIEW_PRIMITIV_LABEL_H
#define BALL_VIEW_PRIMITIV_LABEL_H

#ifndef BALL_VIEW_KERNEL_GEOMETRICOBJECT_H
# include <BALL/VIEW/KERNEL/geometricObject.h>
#endif

#ifndef BALL_VIEW_DATATYPE_VERTEX1_H
#	include <BALL/VIEW/DATATYPE/vertex1.h>
#endif

#include <qfont.h>

namespace BALL
{
	namespace VIEW
	{

		/** Label class.
				An instance of Label represents an instance of the geometric representation "label".
				A label is an information text that can be pinned to another Composite or
				GeometricObject. A label is both visible in the dynamic and static render
				mode of the Scene.
				A label has the following properties. 
				  - color  - the color of the label
					- text   - the text of the label
					- vertex - the position of the label			
				\par
				The class Label is derived from the classes GeometricObject, ColorExtension
				and Vertex. See these classes for further information concerning
				interface and additional methods. \par
				\ingroup ViewPrimitives
		*/
		class BALL_VIEW_EXPORT Label
			: public GeometricObject,
				public Vertex
		{
			public:

			BALL_CREATE(Label)

			/**	@name	Constructors
			*/	
			//@{

			/** Default Constructor.
					Construct new label.
					The properties of this label are set to:
  				  - color  - to the color black
						- text   - to the text "unkown"
		  			- vertex - to the vector (0,0,0)
					\par
					\return      Label new constructed label
					\see         GeometricObject
					\see         Vertex
			*/
			Label()
				throw();

			/** Copy constructor with cloning facility.
					\param       label the label to be copied (cloned)
					\see         GeometricObject
					\see         Vertex
			*/
			Label(const Label& label)
				throw();

			//@}
			/** @name Destructors */
			//@{

			/** Destructor.
			*/
			virtual ~Label()
				throw();

			/** Explicit default initialization.
					Calls GeometricObject::clear
					Calls Vertex::clear
					\see  GeometricObject::clear
					\see  Vertex::clear
			*/
			virtual void clear()
				throw();

			//@}
			/**	@name	Assignment methods */
			//@{

			/** Assignment.
					Assign the label <b> label</b> to this label.
					\param       label the label to be copied
			*/
			void set(const Label& label)
				throw();

			/** Assignment operator.
					Calls set.
			*/
			const Label& operator = (const Label& label)
				throw();

			/** Swapping of label's.
			*/
			void swap(Label& label)
				throw();

			//@}
			/**	@name	Inspectors, Mutators, Accessors */
			//@{

			/** Change the text of the label.
			*/
			void setText(const String& text)
				throw() {text_ = text;}

			/** Inspection of the text of the label.
			*/
			String getText() const
				throw() { return text_;}

			/** Inspection of the expanded text of the label.
			*/
			String getExpandedText() const
				throw();

			///
			const QFont& getFont() const { return font_;}

			///
			void setFont(const QFont& font) { font_ = font;}

			//@}
			/**	@name	debuggers and diagnostics */
			//@{

			/** Internal state and consistency self-validation.
					Initiate self-validation of the internal state and data structure consistencies
					of this label.
					If the internal state of this label is correct (self-validated) and 
					consistent <tt> true</tt> is returned, <tt> false</tt> otherwise. 
					Calls GeometricObject::isValid.
					Calls Vertex::isValid.
					\return			bool <tt> true</tt> if the internal state of this label is correct 
											(self-validated) and consistent, <tt> false</tt> otherwise
					\see        GeometricObject::isValid
					\see        Vertex::isValid
			*/
			virtual bool isValid() const
				throw();

			/** Internal value dump.
					Dump the current value of this label to 
					the output ostream <b> s</b> with dumping depth <b> depth</b>.
					Calls GeometricObject::dump.
					Calls Vertex::dump.
					\param   s output stream where to output the value of this label
					\param   depth the dumping depth
					\see     GeometricObject::dump
					\see     Vertex::dump
			*/
			virtual void dump(std::ostream&  s = std::cout, Size depth = 0) const
				throw();

			protected:
				String text_;
				QFont  font_;

			//@}
		};

	} // namespace VIEW
} // namespace BALL

#endif // BALL_VIEW_PRIMITIV_LABEL_H
