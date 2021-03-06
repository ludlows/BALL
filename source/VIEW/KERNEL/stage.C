#include <BALL/VIEW/KERNEL/stage.h>
#include <BALL/MATHS/plane3.h>
#include <BALL/MATHS/analyticalGeometry.h>

using std::endl;

namespace BALL
{
	namespace VIEW
	{

		LightSource::LightSource()
			throw()
			:	position_(),
				direction_(0, 0, -1),
				angle_(10),
				intensity_(0.8),
				color_(255, 255, 255, 255),
				type_(0),
				relative_(true)	
		{
		}

		LightSource::LightSource(const LightSource& light_source)
			throw()
			: position_(light_source.position_),
				direction_(light_source.direction_),
				angle_(light_source.angle_),
				intensity_(light_source.intensity_),
				color_(light_source.color_),
				type_(light_source.type_),
				relative_(light_source.relative_)
		{
		}

		LightSource& LightSource::operator = (const LightSource& light) 
			throw()
		{
				position_ = light.position_;
			 direction_ = light.direction_;
					 angle_ = light.angle_,
			 intensity_ = light.intensity_;
					 color_ = light.color_,
						type_ = light.type_;
				relative_ = light.relative_;

			return *this;
		}


		bool LightSource::operator == (const LightSource& light_source) const
			throw()
		{
			return position_ 		== light_source.position_ 	&&
						 direction_ 	== light_source.direction_ 	&&
						 angle_ 			== light_source.angle_ 		 	&&
						 intensity_ 	== light_source.intensity_ 	&&
						 color_ 			== light_source.color_ 			&&
						 type_ 				== light_source.type_ 			&&
						 relative_ 		== light_source.relative_;
		}
		
		bool LightSource::operator < (const LightSource& light) const
			throw()
		{
			return this < &light;
		}


		void LightSource::dump(std::ostream& s, Size depth) const
			throw()
		{
			BALL_DUMP_STREAM_PREFIX(s);
			
			BALL_DUMP_DEPTH(s, depth);
			BALL_DUMP_HEADER(s, this, this);

			BALL_DUMP_DEPTH(s, depth);
			s << "Position : " << position_ << std::endl;

			BALL_DUMP_DEPTH(s, depth);
			s << "Direction : " << direction_<< std::endl;

			BALL_DUMP_DEPTH(s, depth);
			s << "Angle : " << angle_ << std::endl;

			BALL_DUMP_DEPTH(s, depth);
			s << "Intensity : " << intensity_ << std::endl;

			BALL_DUMP_DEPTH(s, depth);
			s << "Color : " << color_ << std::endl;

			BALL_DUMP_DEPTH(s, depth);
			s << "Type : " << type_ << std::endl;

			BALL_DUMP_DEPTH(s, depth);
			s << "Relative : " << relative_ << std::endl;

			BALL_DUMP_STREAM_SUFFIX(s);
		}

		Camera::Camera()
			throw()
			: view_point_(0, 0, 0),
				look_at_(0, 0, -1),
				look_up_vector_(0, -1, 0)
		{
			calculateVectors_();
		}

		Camera::Camera(const Camera& camera)
			throw()
			: view_point_(camera.view_point_),
				look_at_(camera.look_at_),
				look_up_vector_(camera.look_up_vector_),
				view_vector_(camera.view_vector_),
				right_vector_(camera.right_vector_)
		{
		}

		Camera::Camera(const Vector3& view_point, const Vector3& look_at, const Vector3& look_up_vector)
			throw()
			: view_point_(view_point),
				look_at_(look_at),
				look_up_vector_(look_up_vector)
		{
			calculateVectors_();
		}

		Camera& Camera::operator = (const Camera& camera) 
			throw()
		{
					 view_point_ = camera.view_point_;
							look_at_ = camera.look_at_;
			 look_up_vector_ = camera.look_up_vector_;
					view_vector_ = camera.view_vector_;
				 right_vector_ = camera.right_vector_;

			return *this;
		}
		
	
		bool Camera::operator < (const Camera& camera) const
			throw()
		{
			return this < &camera;
		}


		bool Camera::operator == (const Camera& camera) const
			throw()
		{
			return 	view_point_ 		== camera.view_point_ &&
							look_at_ 				== camera.look_at_ 		&&
							look_up_vector_ == camera.look_up_vector_;
		}

		void Camera::calculateVectors_()
			throw()
		{
			if (!look_up_vector_.isZero()) look_up_vector_.normalize();
			else look_up_vector_.set(0,1,0);

			view_vector_  = look_at_ - view_point_;
			right_vector_ = look_up_vector_ % view_vector_;
			right_vector_ *= -1.0;

			if (!right_vector_.isZero())	right_vector_.normalize();
			else right_vector_.set(1,0,0);
		}


		void Camera::dump(std::ostream& s, Size depth) const
			throw()
		{
			BALL_DUMP_STREAM_PREFIX(s);
			
			BALL_DUMP_DEPTH(s, depth);
			BALL_DUMP_HEADER(s, this, this);

			BALL_DUMP_DEPTH(s, depth);
			s << "Viewpoint : " << view_point_ << endl;

			BALL_DUMP_DEPTH(s, depth);
			s << "Look at : " << look_at_ << endl;

			BALL_DUMP_DEPTH(s, depth);
			s << "Look up vector: " << look_up_vector_ << endl;

			BALL_DUMP_STREAM_SUFFIX(s);
		}

		String Camera::toString() const
			throw()
		{
			return vector3ToString(view_point_) + " " +
						 vector3ToString(look_at_) + " " +
						 vector3ToString(look_up_vector_);
		}

		bool Camera::readFromString(const String& data)
			throw()
		{
			vector<String> fields;
			if (data.split(fields) != 3) return false;

			Vector3 results[3];
			for (Position p = 0; p < 3; p++)
			{
				if (!stringToVector3(fields[p], results[p])) return false;
			}

			view_point_ = results[0];
			look_at_ 		= results[1];
			look_up_vector_ = results[2];

			calculateVectors_();
			
			return true;
		}
				

		void Camera::rotate(const Quaternion& q, const Vector3& origin)
			throw()
		{
			translate(-origin);
			Matrix4x4  m;
			q.getRotationMatrix(m);

			view_point_     = m*view_point_;
			look_at_        = m*look_at_;
			look_up_vector_ = m*look_up_vector_;

			calculateVectors_();
			translate(origin);
		}

		Stage::Stage()
			throw()
			: background_color_(),
				light_sources_(),
				camera_(),
				show_coordinate_system_(false),
				fog_intensity_(0),
				eye_distance_(2.0),
				focal_distance_(40),
				swap_side_by_side_stereo_(false),
				specular_(0.4),
				diffuse_(0.2),
				ambient_(0.0),
				shininess_(128.0)
		{}

		Stage::Stage(const Stage& stage)
			throw()
			: background_color_(stage.background_color_),
				light_sources_(stage.light_sources_),
				camera_(stage.camera_),
				show_coordinate_system_(false),
				fog_intensity_(stage.fog_intensity_),
				eye_distance_(stage.eye_distance_),
				focal_distance_(stage.focal_distance_),
				swap_side_by_side_stereo_(stage.swap_side_by_side_stereo_),
				specular_(stage.specular_),
				diffuse_(stage.diffuse_),
				ambient_(stage.ambient_),
				shininess_(stage.shininess_)
		{
		}

		void Stage::clear()
			throw()
		{
			background_color_.clear();
			light_sources_.clear();
			camera_ = Camera();
			show_coordinate_system_ = false;
			eye_distance_ = 2.0;
			focal_distance_ = 40;
			swap_side_by_side_stereo_ = false;
			fog_intensity_ = 0;
			specular_ = 0.4;
			diffuse_  = 0.2;
			ambient_  = 0.0;
			shininess_ = 128.0;
		}


		bool Stage::operator == (const Stage& stage) const
			throw()
		{
			return light_sources_ 					== stage.light_sources_ 		&&
						 camera_ 					 				== stage.camera_ 						&&
						 background_color_ 				== stage.background_color_ 	&&
						 show_coordinate_system_ 	== stage.show_coordinate_system_ &&
						 eye_distance_ 						== stage.eye_distance_ 			&&
						 focal_distance_ 					== stage.focal_distance_ 		&&
						 swap_side_by_side_stereo_== stage.swap_side_by_side_stereo_ &&
						 specular_ 								== stage.specular_ &&
						 diffuse_ 								== stage.diffuse_  &&
						 ambient_ 								== stage.ambient_  &&
						 shininess_ 							== stage.shininess_;
		}



		void Stage::dump(std::ostream& s, Size depth) const
			throw()
		{
			BALL_DUMP_STREAM_PREFIX(s);
			
			BALL_DUMP_DEPTH(s, depth);
			BALL_DUMP_HEADER(s, this, this);

			BALL_DUMP_DEPTH(s, depth);
			s << "Light sources: -------------------------------" << endl;

			List<LightSource>::ConstIterator it = light_sources_.begin();
			for (; it != light_sources_.end(); it++)
			{
				it->dump(s, depth + 1);
			}

			BALL_DUMP_DEPTH(s, depth);
			s << "Camera: ---------------------------------------" << endl;
			camera_.dump(s, depth);

			BALL_DUMP_DEPTH(s, depth);
			s << "Background color:  " << background_color_ << endl;

			BALL_DUMP_DEPTH(s, depth);
			s << "Show coordinate system:  " << show_coordinate_system_ << endl;

			BALL_DUMP_DEPTH(s, depth);
			s << "Fog intensity:  " << fog_intensity_ << endl;

			BALL_DUMP_DEPTH(s, depth);
			s << "Eye distance :  " << eye_distance_<< endl;

			BALL_DUMP_DEPTH(s, depth);
			s << "Focal width:  " << focal_distance_<< endl;

			BALL_DUMP_DEPTH(s, depth);
			s << "Swap side by side:  " << swap_side_by_side_stereo_ << endl;

			BALL_DUMP_STREAM_SUFFIX(s);
		}

		Vector3 Stage::calculateRelativeCoordinates(Vector3 pos) const
		{
			// relative in units of right_vector , look_up_vector , view_vector
			// by calculating the normals to three planes
			
			// make sure the three planes are far enough, to be always on one side of them
			const float d = 1000000.0;

			try
			{
				Vector3 dr(camera_.getRightVector() * d);

				Vector3 dv(camera_.getViewVector());
				dv.normalize();
				dv *= d;

				Vector3 du(camera_.getLookUpVector() * d);

				// calculate the planes
				const Plane3 plane_rv(dr, dr);
				const Plane3 plane_uv(du, du);
				const Plane3 plane_vv(dv, dv);

				// distance of the destion of the light source from the three planes
				Vector3 result(
					GetDistance(plane_rv, pos),
					GetDistance(plane_uv, pos),
					GetDistance(plane_vv, pos));
				result -= Vector3(d);

				return -result;
			}
			catch(...)
			{
				Log.error() << "Could not calculate relative light coordinates, degenerated view vectors." << std::endl;
			}

			return Vector3(1.0);
		}

		Vector3 Stage::calculateAbsoluteCoordinates(Vector3 pos) const
		{
			try
			{
				Vector3 dv(camera_.getViewVector());
				dv.normalize();

				return (pos.x * camera_.getRightVector() + 
							  pos.y * camera_.getLookUpVector() +
								pos.z * dv);
			}
			catch(...)
			{
				Log.error() << "Could not calculate absolute light coordinates, degenerated view vectors." << std::endl;
			}

			return Vector3(1.0);
		}

		void Stage::addLightSource(const LightSource& light_source)
			throw()
		{
			light_sources_.push_back(light_source);
		}

		void Stage::removeLightSource(const LightSource& light_source)
			throw()
		{
			List<LightSource>::Iterator it = light_sources_.begin();
			for (; it != light_sources_.end(); it++)
			{
				if (&*it == &light_source) 
				{
					light_sources_.erase(it);
					return;
				}
			}
		}

		void Stage::clearLightSources()
			throw()
		{
			light_sources_.clear();
		}

	} // namespace VIEW
} // namespace BALL
