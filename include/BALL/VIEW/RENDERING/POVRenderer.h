// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// $Id: POVRenderer.h,v 1.8.4.7 2005/09/01 22:18:02 amoll Exp $
//

#ifndef BALL_VIEW_RENDERING_POVRENDERER_H
#define BALL_VIEW_RENDERING_POVRENDERER_H

#ifndef BALL_VIEW_RENDERING_RENDERER_H
# include <BALL/VIEW/RENDERING/renderer.h>
#endif

#ifndef BALL_SYSTEM_FILE_H
# include <BALL/SYSTEM/file.h>
#endif

#ifndef BALL_MATHS_VECTOR3_H
# include <BALL/MATHS/vector3.h>
#endif

#ifndef BALL_MATHS_MATRIX44_H
# include <BALL/MATHS/matrix44.h>
#endif

namespace BALL
{
	namespace VIEW
	{
		class ColorRGBA;
		class ClippingPlane;

		/** POVRenderer class.
		 		This class walks over all the geometric primitives in a Scene
				and exports them into a data file in the POVRay 1.5 format, which can
				be used to render the same scene externally.
				\ingroup ViewRendering
		*/
		class BALL_VIEW_EXPORT POVRenderer : public Renderer
		{
			public:

			BALL_CREATE(POVRenderer)

			struct POVRendererClippingPlane
			{
				public:
					Vector3 normal;
					float translation;
//					 Vector3 translation;
			};

			/** @name Constructors and Destructors.
			 */
			//@{

			/// Default constructor.
			POVRenderer()
				throw();

			/** Detailed constructor.
			 		\param name The name of the file we will create
			 */
			POVRenderer(const String& name)
				throw(Exception::FileNotFound);
			
			// Only for Python
			POVRenderer(const POVRenderer& renderer)
				throw();


			/// Destructor.
			virtual ~POVRenderer()
				throw();

			/// Clear method.
			virtual void clear()
				throw();

			//@}
			/** @name Accessors
			 */
			//@{

			/** Sets the name of the file we will create.
			 		\param name The file name
			 */
			void setFileName(const String& name)
				throw(Exception::FileNotFound);

			/// Set a stream as output device
			void setOstream(std::ostream& out_stream);

			/// 
			void setHumanReadable(bool state)
				throw() { human_readable_ = state;}

			///
			bool isHumanReadable() const
				throw() { return human_readable_;}

			/** Converts a ColorRGBA into a String in POVRay format.
			 */
			String POVColorRGBA(const ColorRGBA& input)
				throw();

			/** Returns the corresponding BALLFinish declaration
			 */
			String POVFinish(const String& object, const ColorRGBA& input)
				throw();

			/** Converts a Vector3 into a String in POVRay format.
			 */
			String POVVector3(Vector3 input)
				throw();

			virtual bool render(const Representation& representation)
				throw();

			//@}
			
			/** @name Processor specific methods
			 */
			//@{
			/** Start method. 
			    This method creates the file and writes the header.
			 */
			virtual bool init(const Stage& stage, float width, float height)
				throw();

			/** Finish method.
			 		This method writes the ending of the file and closes it.
			 */
			virtual bool finish()
				throw();

			void renderSphere_(const Sphere& sphere)
				throw();
			
			void renderDisc_(const Disc& disc)
				throw();

			void renderTube_(const Tube& tube)
				throw();

			void renderTwoColoredTube_(const TwoColoredTube& tube)
				throw();

			void renderMesh_(const Mesh& mesh)
				throw();

			void renderTwoColoredLine_(const TwoColoredLine& line)
				throw();

			void renderLine_(const Line& line)
				throw();

			void renderPoint_(const Point& point)
				throw();

			// do nothing
			void renderLabel_(const Label&)
				throw();

			//@}

			protected:

				std::ostream* outfile_;
				String trimFloatValue_(float value);
				const ColorRGBA& getColor_(const GeometricObject& object);
				void storeColor_(const GeometricObject& object);
				String getColorIndex_(const ColorRGBA& color);

				Vector3   origin_;
				Matrix4x4 rotation_;
				vector<ClippingPlane*> clipping_planes_;
				bool human_readable_;

				typedef HashMap<String, Position> ColorMap;
				ColorMap color_map_;
				vector<const ColorRGBA*> color_vector_;
				vector<const Representation*> representations_;
				HashSet<const Mesh*> wireframes_;
				String font_file_;
				double m_[12];
		};
  
	} // namespace BALL
} // namespace VIEW

#endif // BALL_VIEW_RENDERING_POVRENDERER_H
