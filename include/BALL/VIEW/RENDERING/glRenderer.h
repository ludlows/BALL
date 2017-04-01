// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// $Id: glRenderer.h,v 1.32.2.13 2005/12/07 15:50:04 amoll Exp $
//

#ifndef BALL_VIEW_RENDERING_GLRENDERER_H
#define BALL_VIEW_RENDERING_GLRENDERER_H

#ifndef BALL_VIEW_RENDERING_RENDERER_H
#	include <BALL/VIEW/RENDERING/renderer.h>
#endif

#ifndef BALL_MATHS_QUATERNION_H
# include <BALL/MATHS/quaternion.h>
#endif 

#ifndef BALL_VIEW_DATATYPE_COLORRGBA_H
# include <BALL/VIEW/DATATYPE/colorRGBA.h>
#endif

#ifndef BALL_VIEW_KERNEL_GEOMETRICOBJECT_H
#	include <BALL/VIEW/KERNEL/geometricObject.h>
#endif

#ifndef BALL_VIEW_KERNEL_STAGE_H
# include <BALL/VIEW/KERNEL/stage.h>
#endif

#ifndef BALL_VIEW_RENDERING_GLDISPLAYLIST_H
# include <BALL/VIEW/RENDERING/glDisplayList.h>
#endif

class QFont;

namespace BALL
{
// defines the maximal number of GL-objects, which can be selected in picking mode
// a number as big as 100.000 is needed for large molecules, just to be sure we use a million
#define BALL_GLRENDERER_PICKING_NUMBER_OF_MAX_OBJECTS 1000000
	namespace VIEW
	{
		class GLDisplayList;
		class Scene;
		class MeshBuffer;

		/** GLRenderer
		 		Renderer which provides hardware accelerated OPENGL rendering.
				\ingroup ViewRendering
		*/
		class BALL_VIEW_EXPORT GLRenderer
			: public Renderer
		{
			friend class Scene;
			public:

			///
 			enum StereoMode
			{
				NO_STEREO = 0,

				/// Stereo mode for shutter glasses
				ACTIVE_STEREO,

				/// Stereo mode for output on two projectors
				DUAL_VIEW_STEREO
			};

			///
			enum RenderMode
			{
				///
				RENDER_MODE_UNDEFINED = 0,
				
				///
				RENDER_MODE_SOLID,

				///
				RENDER_MODE_TRANSPARENT,

				///
				RENDER_MODE_ALWAYS_FRONT
			};

			/** WrongModes Exception class.
					This exeption will be thrown if the <b> drawing_precision_</b> or
					<b> drawing_mode_</b> are not allowed. 
					\see         GeneralException			
			*/
			class BALL_VIEW_EXPORT WrongModes:	public Exception::GeneralException
			{
				public:

				WrongModes(const char* file, int line, int mode, int precision)
					throw();
			};

			/** @name Type Definitions
			*/
			//@{
			
			/// Typedef for OPENGL names
			typedef unsigned int Name;
			
			//@}
			/**	@name	Constructors and Destructors
			*/	
			//@{

			/** Default Constructor.
			*/
			GLRenderer()
				throw();

			/** Destructor
			*/
			virtual ~GLRenderer()
				throw();

			/** Explicit default initialization.
			*/
			virtual void clear()
				throw();

			//@}
			/**	@name	Accessors: inspectors and mutators 
			*/
			//@{		

			///
			void dump(std::ostream& s, Size depth) const
				throw();

			///
			inline Name getName(const GeometricObject& object)
				throw();

			///
			GeometricObject* getObject(GLRenderer::Name name) const
				throw();

			/** Initialise the renderer, by calling the init method below
			 		This method is called by Scene::initializeGL.
			*/
			virtual bool init(Scene& scene)
				throw();

			/// Initialise the renderer, e.g. the display lists.
			virtual bool init(const Stage& stage, float height, float width)
				throw();

			/// Set the light sources according to the stage
			virtual void setLights(bool reset_all = false)
				throw();

			/** Pick geometric objects
			 		\param x1, y1, x2, y2 the rectangle of the selection
			*/
			void pickObjects1(Position x1, Position y1, Position x2, Position y2)
				throw();

			/** Pick geometric objects method2.
			 		Call this method after pickObjects1 and rendering the representations.
					\param objects returns the picked objects
			*/
			void pickObjects2(List<GeometricObject*>& objects)
				throw();

			///
			void enterPickingMode();

			///
			void exitPickingMode();

			/**
			 */
			void setSize(float width, float height)
				throw();

			///
			float getXScale() const
				throw();

			///
			float getYScale() const
				throw();

			/** Update the camera position with gluLookAt,
			 		either from a given Camera, or from the default
					Stage.
			*/
			void updateCamera(const Camera* camera = 0)
				throw();

			/// Update the background color from the stage
			void updateBackgroundColor()
				throw();

			// Initialise transparent rendering
			void initTransparent() 
				throw();

			// Initialise solid rendering
			void initSolid()
				throw();
			
			// Initialise always front rendering
			void initAlwaysFront()
				throw();
			
			/// Remove all VertexBuffer and DisplayLists for the given Representation
			void removeRepresentation(const Representation& rep)
				throw();

			/// Buffer the visualisation for the given Representation into OpenGL VertexBuffer Objects and DisplayLists.
			void bufferRepresentation(const Representation& rep)
				throw();

			/// Draw the visualisation of the given Representation from the VertexBuffers and a DisplayList.
			void drawBuffered(const Representation& rep)
				throw();

			/// Test if a Representation has a DisplayList.
			bool hasDisplayListFor(const Representation& rep) const
				throw();
			
			///
			void setStereoMode(StereoMode state)
				throw();

			///
			StereoMode getStereoMode() const
				throw();

			///
			RenderMode getRenderMode() const
				throw();

			///
			void setRenderMode(RenderMode mode) { render_mode_ = mode;}
			
			///
			virtual bool render(const Representation& representation, bool for_display_list = false)
				throw();

			/** Test if a given opengl extension is supported by the current driver.
			 		Call this only after Scene::initializeGL();
			*/
			bool isExtensionSupported(const String& extension) const
				throw();

			/// 
			void clearVertexBuffersFor(Representation& rep)
				throw();

			///
			bool vertexBuffersSupported() const
				throw();

			///
			String getVendor();

			///
			String getRenderer();

			///
			String getOpenGLVersion();

			///
			vector<String> getExtensions();

			///
			bool enableVertexBuffers(bool state)
				throw();

			///
			bool vertexBuffersEnabled() const;

			///
			DrawingMode getDrawingMode() const;

			///
			void initPerspective();

			//@}
			protected:

			///
			virtual void renderLabel_(const Label& /*label*/)
				throw();

			///
			virtual void renderLine_(const Line& /*line*/)
				throw();

			///
			virtual void renderMesh_(const Mesh& /*mesh*/)
				throw();

			///
			void initDrawingMeshes_();

			///
			void initDrawingOthers_();

			///
			virtual void renderPoint_(const Point& /*point*/)
				throw();

			///
			virtual void renderSimpleBox_(const SimpleBox& /*box*/)
				throw();

			///
			virtual void renderBox_(const Box& /*box*/)
				throw();

			///
			virtual void renderSphere_(const Sphere& /*sphere*/)
				throw();

			///
			virtual void renderDisc_(const Disc& /*disc*/)
				throw();

			///
			virtual void renderTube_(const Tube& /*tube*/)
				throw();

			///
			virtual void renderTwoColoredLine_(const TwoColoredLine& /*two_colored_line*/)
				throw();

			///
			virtual void renderTwoColoredTube_(const TwoColoredTube& /*two_colored_tube*/)
				throw();

			///
			virtual void renderClippingPlane_(const ClippingPlane& plane)
				throw();

			//_
			void setColor4ub_(const GeometricObject& object)
				throw();

			//_
			void createSpheres_()
				throw();
			
			//_
			void createTubes_()
				throw();

			//_
			void createBoxes_()
				throw();

			//_
			void createDottedSphere_(int precision)
				throw();
			
			//_
			void subdivideTriangle_(Vector3& v1, Vector3& v2, Vector3& v3, int precision)
				throw();

			//_
			void createLineBox_()
				throw();

			//_
			void createDotBox_()
				throw();

			//_
			void createSolidBox_()
				throw();

			//_
			void clearNames_()
				throw();

			//_
			void normalVector3_(const Vector3& v) 
				throw();

			//_
			void vertexVector3_(const Vector3& v)
				throw();

			//_
			void translateVector3_(const Vector3& v)
				throw();

			//_
			void scaleVector3_(const Vector3& v)
				throw();

			//_
			void rotateVector3Angle_(const Vector3& v, Real angle)
				throw();

			//_
			void scale_(float f)
				throw();

			//_
			void setColorRGBA_(const ColorRGBA& color)
				throw();

			void initGLU_(DrawingMode mode);

			//_
 			GLubyte* generateBitmapFromText_(const String& text, const QFont& font, int& width, int& height) const
				throw();

			Scene* 								scene_;

			///
			DrawingMode 					drawing_mode_;

			///
			Index 								drawing_precision_;

			//_
			float 								x_scale_;

			//_
			float 								y_scale_;

			GLDisplayList* 				GL_spheres_list_;
			GLDisplayList* 				GL_tubes_list_;
			GLDisplayList* 				GL_boxes_list_;
			GLDisplayList* 				sphere_list_;

			/* static array of vertices for sphere dots */
			static const float sphere_vertices_[12][3];
			static const int 		sphere_indices_[20][3];

			// naming of geometric objects
			typedef HashMap<const GeometricObject*, Name> NameHashMap;
			typedef HashMap<Name, const GeometricObject*> GeometricObjectHashMap;
			typedef HashMap<const Representation*, GLDisplayList*> DisplayListHashMap;
			typedef HashMap<const Representation*, vector<MeshBuffer*> > MeshBufferHashMap;

			GeometricObjectHashMap	name_to_object_;
			NameHashMap							object_to_name_;
			DisplayListHashMap 			display_lists_;
			MeshBufferHashMap 			rep_to_buffers_;
			Name 										all_names_;
			GLuint 									object_buffer_[BALL_GLRENDERER_PICKING_NUMBER_OF_MAX_OBJECTS];
			Vector3 								normal_vector_;
			ColorRGBA 							dummy_color_;
			const ColorRGBA* 				last_color_;

			StereoMode 							stereo_;
			RenderMode 							render_mode_;

			bool 										use_vertex_buffer_;
			bool 										picking_mode_;
			ModelType 							model_type_;
			Position 								display_lists_index_;
			bool 										single_pick_;
			bool 										drawed_other_object_;
			bool 										drawed_mesh_;
			GLUquadricObj*  GLU_quadric_obj_;
		};

#	ifndef BALL_NO_INLINE_FUNCTIONS
#		include <BALL/VIEW/RENDERING/glRenderer.iC>
#	endif

	} // namespace VIEW
} // namespace BALL

#endif // BALL_VIEW_RENDERING_GLRENDERER_H
