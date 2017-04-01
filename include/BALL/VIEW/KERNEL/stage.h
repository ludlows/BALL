// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// $Id: stage.h,v 1.18.4.5 2005/11/03 17:37:22 amoll Exp $

#ifndef BALL_VIEW_KERNEL_STAGE_H
#define BALL_VIEW_KERNEL_STAGE_H

#ifndef BALL_MATHS_VECTOR3_H
# include <BALL/MATHS/vector3.h>
#endif

#ifndef BALL_MATHS_QUATERNION_H
#	include <BALL/MATHS/quaternion.h>
#endif

#ifndef BALL_MATHS_MATRIX44_H
#	include <BALL/MATHS/matrix44.h>
#endif

#ifndef BALL_MATHS_ANGLE_H
# include <BALL/MATHS/angle.h>
#endif

#ifndef BALL_VIEW_DATATYPE_COLORRGBA_H
# include <BALL/VIEW/DATATYPE/colorRGBA.h>
#endif

#ifndef BALL_DATATYPE_LIST_H
# include <BALL/DATATYPE/list.h>
#endif

#ifndef  BALL_VIEW_KERNEL_REPRESENTATION_H
# include <BALL/VIEW/KERNEL/representation.h>
#endif

namespace BALL
{
	namespace VIEW
	{
		/** Light source is mainly used for Renderer classes (e.g. OpenGL and POVRay).
		 		Currently we support ambient, positional and directional light sources.
		 		The Position and direction of lights can be stored twofold:
				- Either with absolute room coordinates and a direction vector or \\
				- Relative to a Camera.
				In this case, the position and direction vector are stored as multiples of
				look right, look up and view vector.
				\ingroup  ViewKernelStage
		 */
		class BALL_VIEW_EXPORT LightSource
		{
			public:

			/// Enumeration of different types of lights
			enum Types
			{

			 	/**	Ambient light doesn't come from any particular direction. 
						All the objects in the scene will be lit up by the ambient light.
				*/
				AMBIENT = 0,

				/** Diffuse light is created the light source and is reflected off the surface 
						of any object in the scene. Any surface of an object that the light hits 
						directly will be very bright, and areas the light barely gets to will be darker.
				*/
				POSITIONAL,

				/**
				*/
				DIRECTIONAL
			};

			/**	@name	Constructors and Destructors
			*/	
			//@{

			/** Constructor
			 */
			LightSource()
				throw();

			/** Copy Constructor
			 */
			LightSource(const LightSource& light_source)
				throw();

			/** Destructor
			 */
			virtual ~LightSource()
				throw() {}

			//@}
			/**	@name	Accessors
			*/
			//@{

			/// Get position
			const Vector3& getPosition() const
				throw() { return position_;}

			/// Set position
			void setPosition(const Vector3& position)
				throw() { position_ = position; }

			/// Get the direction vector of the light
			const Vector3& getDirection() const
				throw() { return direction_;}

			/// Set the direction vector of the light
			void setDirection(const Vector3& direction)
				throw() { direction_ = direction;}

			/// Get the angle of the light cone
			const Angle& getAngle() const
				throw() { return angle_;}

			/// Set the angle of the light cone
			void setAngle(const Angle& angle)
				throw() { angle_ = angle;}
			
			/** Get the light intensity.
					0 is the minumum, 1 is the maximum.
			*/
			float getIntensity() const
				throw() { return intensity_;}

			/** Set the intensity.
					0 is the minumum, 1 is the maximum.
			*/
			void setIntensity(float intensity)
				throw() { intensity_ = intensity;}
				
			/** Get the color of the light.
			 		The alpha channel of the color is ignored.
			*/
			const ColorRGBA& getColor() const
				throw() { return color_;}

			/** Set the color of the light.
			 		The alpha channel of the color is ignored.
			*/
			void setColor(const ColorRGBA& color)
				throw() { color_ = color;}
			
			/** Get the type of the light.
			 		\see Types
			*/
			Index getType() const
				throw() { return type_;}

			/** Set the type of the light.
			 		\see Types
			*/
			void setType(Types type)
				throw() { type_ = type;}

			/// If set to true, the LightSource will move with the Camera
			void setRelativeToCamera(bool state)
				throw() { relative_ = state;}

			/// Test if a LightSource will move with the Camera
			bool isRelativeToCamera() const
				throw() { return relative_;}

			///
			LightSource& operator = (const LightSource& light) throw();

			/// needed for MSVC, dont use it otherwise!
			bool operator < (const LightSource& light) const
				throw();

			//@}
			/**	@name Predicates
			*/
			//@{
			
			///
			bool operator == (const LightSource& light_source) const
				throw();

			//@}
			
			/** Internal value dump.
					Dump the current state of this instance to 
					the output ostream <b> s</b> with dumping depth <b> depth</b>.
					\param   s output stream 
					\param   depth the dumping depth
			*/
			virtual void dump(std::ostream& s = std::cout, Size depth = 0) const
				throw();

			protected:

			//_ Position of the light source
			Vector3 		position_;

			//_ Direction of the light cone
			Vector3 		direction_;

			//_
			Vector3     r_position_;

			//_
			Vector3  		r_direction_;

			//_ Angle of the light cone
			Angle 			angle_;

			//_ Intensity of the light
			float 			intensity_;

			//_ Color of the light
			ColorRGBA 	color_;

			//_ Enumeration of types for further usage
			Index 			type_;

			//_ Relative to camera
			bool 				relative_;
		};
				

		/** Camera with viewpoint, a look at point and an up-vector.
				\ingroup ViewKernelStage
		 */
		class BALL_VIEW_EXPORT Camera
		{
			public:

			/**	@name	Constructors and Destructors
			*/	
			//@{
			
			/// Constructor
			Camera()
				throw();

			/// Copy Constructor
			Camera(const Camera& camera)
				throw();

			Camera(const Vector3& view_point, const Vector3& look_at, const Vector3& look_up_vector)
				throw();

			/// Destructor
			virtual ~Camera()
				throw() {}

			///
			Camera& operator = (const Camera& camera) throw();

			//@}
			/**	@name	Accessors
			*/
			//@{

			/// Get the position of the camera
			const Vector3& getViewPoint() const
				throw() { return view_point_;}

			/// Set the position of the camera
			void setViewPoint(const Vector3& view_point) 
				throw() { view_point_ = view_point; calculateVectors_();}

			/// Get the direction of the camera
			const Vector3& getLookAtPosition() const
				throw() { return look_at_;}

			/// Set the direction of the camera
			void setLookAtPosition(const Vector3& look_at) 
				throw() { look_at_ = look_at; calculateVectors_();}

			/// Get the look up vector
			const Vector3& getLookUpVector() const
				throw() { return look_up_vector_;}

			/// Set the look up vector
			void setLookUpVector(const Vector3& look_up_vector) 
				throw() { look_up_vector_ = look_up_vector; calculateVectors_();}

			/// Get the distance between the view point and the look at point
			float getDistance() const
				throw() { return view_point_.getDistance(look_at_);}

			/// Get the view vector
			Vector3 getViewVector() const
				throw() { return view_vector_;}

			/// Get an vector orthogonal to the viewing vector and showing to the right
			Vector3 getRightVector() const
				throw() { return right_vector_;}
			
			/// Translate the view point and the point the camera is looking to by a given vector
			void translate(const Vector3& v) 
				throw() { view_point_ += v; look_at_ += v; calculateVectors_();}

			/// Rotate the camera
			void rotate(const Quaternion& q, const Vector3& origin)
				throw();

			/// Reset Camera to standard values
			virtual void clear()
				throw() { *this = Camera();}

			///
			String toString() const
				throw();

			///
			bool readFromString(const String& data)
				throw();

			//@}
			/**	@name Predicates
			*/
			//@{
			
			///
			bool operator == (const Camera& camera) const
				throw();

			/// Needed for MSVC
			bool operator < (const Camera& camera) const
				throw();

			//@}
	
			/** Internal value dump.
					Dump the current state of this instance to 
					the output ostream <b> s</b> with dumping depth <b> depth</b>.
					\param   s output stream 
					\param   depth the dumping depth
			*/
			virtual void dump(std::ostream& s = std::cout, Size depth = 0) const
				throw();

			protected:

			//_
			void calculateVectors_()
				throw();

			//_
			Vector3 					view_point_;

			//_
			Vector3 					look_at_;

			//_
			Vector3 					look_up_vector_;

			/*_ The viewing vector.
			 		Stored only for better performance in the scene.
			*/
			Vector3 					view_vector_;
	
			/*_ The orthogonal vector to the viewing vector, which is showing to the right.
			 		Stored only for better performance in the scene.
			*/
			Vector3 					right_vector_;
		};

		/** A Stage has a Camera, LightSources and a background color.
		 		It stores also the eye distance for the stereo view.
		 		Finally a flag can be set, so that a coordinate system will be shown.
				\ingroup ViewKernelStage
		*/
		class BALL_VIEW_EXPORT Stage
		{
			public:

			/**	@name	Constructors and Destructors
			*/	
			//@{

			/** Default Constructor
			*/
			Stage()
				throw();
			
			///	Copy constructor
			Stage(const Stage& stage)
				throw();

			/// Destructor
			virtual ~Stage()
				throw() {}

			/// Explicit default initialization.
			virtual void clear()
				throw();

			//@}
			/**	@name	Accessors
			*/
			//@{
			
			/// Get the light sources (const)
			virtual const List<LightSource>& getLightSources() const
				throw() { return light_sources_;}

			/// Add a light source
			virtual void addLightSource(const LightSource& light_source)
				throw();

			/// Remove a light source
			virtual void removeLightSource(const LightSource& light_source) 
				throw();

			///
			void clearLightSources()
				throw();
			
			/// Get the camera
			virtual Camera& getCamera() 
				throw() { return camera_;}

			/// Get the camera (const)
			virtual const Camera& getCamera() const
				throw() { return camera_;}

			/** Set the camera of the stage
			 */
			virtual void setCamera(const Camera& camera)
				throw() { camera_ = camera;}

			/// Get the background color
			virtual const ColorRGBA& getBackgroundColor() const
				throw() { return background_color_;}

			/// Set the background color
			virtual void setBackgroundColor(const ColorRGBA& color)
				throw() { background_color_ = color;}

			/// Show coordinate system
			void showCoordinateSystem(bool state)
				throw() { show_coordinate_system_ = state;}

			/// Shows coordinate system
			bool coordinateSystemEnabled() const
				throw() { return show_coordinate_system_;}

			/// Set the eye distance for the stereo view
			void setEyeDistance(float value) 
				throw() { eye_distance_ = value;}

			/// Get the eye distance for the stereo view
			float getEyeDistance() const
				throw() { return eye_distance_;}
				
			/// Set the focal distance for the stereo view
			void setFocalDistance(float value) 
				throw() { focal_distance_ = value;}

			/// Get the focal distance for the stereo view
			float getFocalDistance() const
				throw() { return focal_distance_;}

			/// Settings for side by side stereo side swapping
			void setSwapSideBySideStereo(bool state)
				throw() { swap_side_by_side_stereo_ = state;}

			/// Get settings for side by side stereo side swapping
			bool swapSideBySideStereo() const
				throw() { return swap_side_by_side_stereo_;}

			///
			float getFogIntensity() const
				throw() { return fog_intensity_;}

			///
			void setFogIntensity(float value)
				throw() { fog_intensity_ = value;}
				
			///
			float getSpecularIntensity() const
				throw() { return specular_;}

			///
			void setSpecularIntensity(float value)
				throw() { specular_ = value;}
					
			///
			float getDiffuseIntensity() const
				throw() { return diffuse_;}

			///
			void setDiffuseIntensity(float value)
				throw() { diffuse_ = value;}
					
			///
			float getAmbientIntensity() const
				throw() { return ambient_;}

			///
			void setAmbientIntensity(float value)
				throw() { ambient_ = value;}
					
			///
			float getShininess() const
				throw() { return shininess_;}

			///
			void setShininess(float value)
				throw() { shininess_ = value;}
			
			//@}
			/**	@name Predicates
			*/
			//@{
			
			///
			bool operator == (const Stage& stage) const
				throw();

			/** Calculate coordiantes relative to the position of the camera in units of 
			 		right_vector, look_up_vector and view_vector.
					This is done by calculating the normals to three planes, spaned by these three vectors.
					This method is e.g. used to store the coordinates of the relative light sources in the INIFile,
					or in the LightSettings dialog.
					@return Vector3(times right_vector, times look_up_vector, times view_vector)
			*/
			Vector3 calculateRelativeCoordinates(Vector3 pos) const;

			/** Calculate absolute room coordinates from relative coordinates.
			 		@see calculateRelativeCoordinates
			*/
			Vector3 calculateAbsoluteCoordinates(Vector3 pos) const;

			//@}
			
			/** Internal value dump.
					Dump the current state of this instance to 
					the output ostream <b> s</b> with dumping depth <b> depth</b>.
					\param   s output stream 
					\param   depth the dumping depth
			*/
			virtual void dump(std::ostream& s = std::cout, Size depth = 0) const
				throw();

			protected:

			//_
			ColorRGBA 					background_color_;

			//_
			List<LightSource> 	light_sources_;

			//_
			Camera 						 	camera_;

			//_
			bool 								show_coordinate_system_;

			//_
			float 							fog_intensity_;

			//_
			float 							eye_distance_;
			
			//_
			float 							focal_distance_;

			//_
			bool 								swap_side_by_side_stereo_;

			float 							specular_;
			float 							diffuse_;
			float 							ambient_;
			float 							shininess_;
		};

	} // namespace VIEW
} // namespace BALL

#endif // BALL_VIEW_KERNEL_STAGE_H
