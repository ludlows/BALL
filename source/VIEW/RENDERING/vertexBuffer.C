// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// $Id: vertexBuffer.C,v 1.5.2.6 2005/11/22 16:28:00 amoll Exp $

#include <BALL/CONFIG/config.h>
// prevent typedef clash under Linux
#ifdef BALL_ENABLE_VERTEX_BUFFER

#define QT_CLEAN_NAMESPACE
#include <qgl.h>

#ifdef _WINDOWS
// Header Files For Windows
 #define WINDOWS_LEAN_AND_MEAN
 #include <windows.h>
 #include <wingdi.h>	
#else
 #define GLX_GLXEXT_PROTOTYPES // required for Mesa-like implementations
 #include <GL/gl.h>
 #include <GL/glx.h>
 #include <GL/glext.h>
#endif

#include <BALL/VIEW/RENDERING/vertexBuffer.h>
#include <BALL/VIEW/RENDERING/glRenderer.h>
#include <BALL/VIEW/PRIMITIVES/mesh.h>
#include <BALL/VIEW/KERNEL/common.h>

#ifndef APIENTRY
#define APIENTRY
#endif
#ifndef APIENTRYP
#define APIENTRYP APIENTRY *
#endif

#ifndef GLAPI
#define GLAPI extern
#endif

namespace BALL
{
	namespace VIEW
	{

	// declare gl methods pointer
	GLAPI void glBindBuffer (GLenum, GLuint);
	GLAPI void glDeleteBuffers (GLsizei, const GLuint *);
	GLAPI void glGenBuffers (GLsizei, GLuint *);
	GLAPI void glBufferData (GLenum, GLsizei*, const GLvoid *, GLenum);
	GLAPI GLvoid* APIENTRY glMapBuffer (GLenum, GLenum);
	GLAPI GLboolean APIENTRY glUnmapBuffer (GLenum);

	// VBO Extension Definitions, From glext.h
 #define GL_ARRAY_BUFFER_ARB 0x8892
 #define GL_STATIC_DRAW_ARB 0x88E4
 #define GL_ELEMENT_ARRAY_BUFFER_ARB  0x8893

	typedef void (APIENTRY * PFNGLBINDBUFFERARBPROC) (GLenum target, GLuint buffer);
	typedef void (APIENTRY * PFNGLDELETEBUFFERSARBPROC) (GLsizei n, const GLuint *buffers);
	typedef void (APIENTRY * PFNGLGENBUFFERSARBPROC) (GLsizei n, GLuint *buffers);
	typedef void (APIENTRY * PFNGLBUFFERDATAARBPROC) (GLenum target, int size, const GLvoid *data, GLenum usage);

	// VBO Extension Function Pointers
	PFNGLGENBUFFERSARBPROC glGenBuffersARB = NULL;					// VBO Name Generation Procedure
	PFNGLBINDBUFFERARBPROC glBindBufferARB = NULL;					// VBO Bind Procedure
	PFNGLBUFFERDATAARBPROC glBufferDataARB = NULL;					// VBO Data Loading Procedure
	PFNGLDELETEBUFFERSARBPROC glDeleteBuffersARB = NULL;			// VBO Deletion Procedure

		bool MeshBuffer::initGL()
		{
			if (glGenBuffersARB != 0) return true;
			// obtain gl method pointers
		#ifdef _WINDOWS
				glGenBuffersARB = (PFNGLGENBUFFERSARBPROC) wglGetProcAddress("glGenBuffersARB");
				glBindBufferARB = (PFNGLBINDBUFFERARBPROC) wglGetProcAddress("glBindBufferARB");
				glBufferDataARB = (PFNGLBUFFERDATAARBPROC) wglGetProcAddress("glBufferDataARB");
				glDeleteBuffersARB = (PFNGLDELETEBUFFERSARBPROC) wglGetProcAddress("glDeleteBuffersARB");
		#else
			#ifdef GLX_ARB_get_proc_address
				glGenBuffersARB = (PFNGLGENBUFFERSARBPROC) glXGetProcAddressARB((const GLubyte*)"glGenBuffersARB");
				glBindBufferARB = (PFNGLBINDBUFFERARBPROC) glXGetProcAddressARB((const GLubyte*)"glBindBufferARB");
				glBufferDataARB = (PFNGLBUFFERDATAARBPROC) glXGetProcAddressARB((const GLubyte*)"glBufferDataARB");
				glDeleteBuffersARB = (PFNGLDELETEBUFFERSARBPROC) glXGetProcAddressARB((const GLubyte*)"glDeleteBuffersARB");
			#endif
		#endif

			return (glGenBuffersARB != 0);
		}


		GLRenderer* MeshBuffer::gl_renderer_ = 0;

		MeshBuffer::MeshBuffer()
		: mesh_(0),
			buffer_(),
			filled_(false),
			busy_(false),
			multiple_colors_(false),
			color_("0000FF"),
			vertices_(0),
			triangles_(0)
		{
			buffer_[0] = buffer_[1] = buffer_[2] = buffer_[3] = 0;
		}

		MeshBuffer::MeshBuffer(const MeshBuffer& mesh_buffer)
		: mesh_(mesh_buffer.mesh_),
			buffer_(),
			filled_(false),
			busy_(false),
			multiple_colors_(false),
			color_("0000FF")
		{
			buffer_[0] = buffer_[1] = buffer_[2] = buffer_[3] = 0;
		}

		const MeshBuffer& MeshBuffer::operator = (const MeshBuffer& mesh_buffer)
			throw()
		{
			mesh_ = mesh_buffer.mesh_;
			filled_ = false;
			buffer_[0] = buffer_[1] = buffer_[2] = buffer_[3] = 0;
			return *this;
		}

		bool MeshBuffer::initialize()
		{
			if (mesh_ == 0 || busy_) return false;

			busy_ = true;

			clearBuffer();

			vertices_  = mesh_->vertex.size();
			triangles_ = mesh_->triangle.size();
			
			// colors, normals, indices, vertex
			// Get valid Names
			glGenBuffersARB(4, buffer_);

			////////////////////////////////////////////////////////////
			/// upload colors, normals and vertices
			////////////////////////////////////////////////////////////
			float* data = new float[vertices_ * 4];
			if (data == 0) return false;
			for (Size index = 0; index < vertices_; ++index)
			{
				const Size start = index * 3;
				data[start] = mesh_->vertex[index].x;
				data[start + 1] = mesh_->vertex[index].y;
				data[start + 2] = mesh_->vertex[index].z;
			}

			glBindBufferARB(GL_ARRAY_BUFFER_ARB, buffer_[0]);
			glBufferDataARB(GL_ARRAY_BUFFER_ARB, sizeof(float) * vertices_ * 3, data, GL_STATIC_DRAW_ARB);

			for (Size index = 0; index < vertices_; ++index)
			{
				const Size start = index * 3;
				data[start] = mesh_->normal[index].x;
				data[start + 1] = mesh_->normal[index].y;
				data[start + 2] = mesh_->normal[index].z;
			}

			glBindBufferARB(GL_ARRAY_BUFFER_ARB, buffer_[1]);
			glBufferDataARB(GL_ARRAY_BUFFER_ARB, sizeof(float) * vertices_ * 3, data, GL_STATIC_DRAW_ARB);

			if (mesh_->colors.size() > 1)
			{
				multiple_colors_ = true;
				for (Size index = 0; index < vertices_; ++index)
				{
					const Size start = index * 4;
					data[start] = (float) mesh_->colors[index].getRed();
					data[start + 1] = (float) mesh_->colors[index].getGreen();
					data[start + 2] = (float) mesh_->colors[index].getBlue();
					data[start + 3] = (float) mesh_->colors[index].getAlpha();
				}

				glBindBufferARB(GL_ARRAY_BUFFER_ARB, buffer_[2]);
				glBufferDataARB(GL_ARRAY_BUFFER_ARB, sizeof(float) * vertices_ * 4, data, GL_STATIC_DRAW_ARB);
			}
			else
			{
				multiple_colors_ = false;
				if (mesh_->colors.size() == 1)
				{
					color_ = mesh_->colors[0];
				}
			}

			delete[] data;

			////////////////////////////////////////////////////////////
			/// upload triangle indices
			////////////////////////////////////////////////////////////
			unsigned int* indices = new unsigned int[triangles_ * 3];
			if (indices == 0) return false;
			for (Size index = 0; index < triangles_; ++index)
			{
				const Size start = index * 3;
				indices[start] = mesh_->triangle[index].v1;
				indices[start + 1] = mesh_->triangle[index].v2;
				indices[start + 2] = mesh_->triangle[index].v3;
			}

			glBindBufferARB(GL_ELEMENT_ARRAY_BUFFER_ARB, buffer_[3]);
			glBufferDataARB(GL_ELEMENT_ARRAY_BUFFER_ARB, sizeof(unsigned int) * triangles_ * 3, 
											indices, GL_STATIC_DRAW_ARB);
			delete[] indices;

			////////////////////////////////////////////////////////////
			/// cleanup
			////////////////////////////////////////////////////////////
			glBindBufferARB(GL_ARRAY_BUFFER_ARB, 0);
			glBindBufferARB(GL_ELEMENT_ARRAY_BUFFER_ARB, 0);

			filled_ = true;
			busy_ = false;
			return true;
		}

		void MeshBuffer::clear()
		{
			while (busy_)
			{
		#ifdef BALL_PLATFORM_WINDOWS
				Sleep(1);
		#else
				// sleep(1);
		#endif
			}
			mesh_ = 0;
			clearBuffer();
		}

		void MeshBuffer::clearBuffer()
		{
			if (!filled_) return;

			glDeleteBuffersARB(4, buffer_);
			filled_ = false;
		}

		bool MeshBuffer::draw()
		{
			if (!filled_ || gl_renderer_ == 0 || busy_) return false;

			busy_ = true;

			glEnableClientState(GL_VERTEX_ARRAY);
			glEnableClientState(GL_INDEX_ARRAY);

			// vertices
			glBindBufferARB(GL_ARRAY_BUFFER_ARB, buffer_[0]);
			glVertexPointer(3, GL_FLOAT, 0, 0); 

			// indices
			glBindBufferARB(GL_ELEMENT_ARRAY_BUFFER_ARB, buffer_[3]);
			glIndexPointer(GL_UNSIGNED_INT, 0, 0);

			// colors
			if (multiple_colors_)
			{
				glEnableClientState(GL_COLOR_ARRAY);
				glBindBufferARB(GL_ARRAY_BUFFER_ARB, buffer_[2]);
				glColorPointer (4, GL_FLOAT, 0, 0);
			}
			else
			{
				glColor4ub(color_.getRed(), color_.getGreen(), color_.getBlue(), color_.getAlpha());
			}


			////////////////////////////////////////////////////////////
			/// draw it
			////////////////////////////////////////////////////////////
			DrawingMode drawing_mode = gl_renderer_->getDrawingMode();
			if (drawing_mode == DRAWING_MODE_SOLID)
			{
				glEnableClientState(GL_NORMAL_ARRAY);
				glBindBufferARB(GL_ARRAY_BUFFER_ARB, buffer_[1]);
				glNormalPointer(GL_FLOAT, 0, 0);

				glDrawElements(GL_TRIANGLES, triangles_ * 3, GL_UNSIGNED_INT, 0);
				glDisableClientState(GL_NORMAL_ARRAY);
			}
			else if (drawing_mode == DRAWING_MODE_WIREFRAME)
			{
				// do nothing
			}
			else // dots
			{
				glDrawElements(GL_POINTS, vertices_ * 3, GL_UNSIGNED_INT, 0);
			}

			////////////////////////////////////////////////////////////
			/// cleanup
			////////////////////////////////////////////////////////////
			glBindBufferARB(GL_ARRAY_BUFFER_ARB, 0);
			glBindBufferARB(GL_ELEMENT_ARRAY_BUFFER_ARB, 0);
			glDisableClientState(GL_VERTEX_ARRAY);
			glDisableClientState(GL_COLOR_ARRAY);
			glDisableClientState(GL_INDEX_ARRAY);
			
			busy_ = false;
			return true;
		}

		MeshBuffer::~MeshBuffer()
		{
			clear();
		}

	} 

}

#endif
