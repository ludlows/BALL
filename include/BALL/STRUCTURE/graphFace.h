// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// $Id: graphFace.h,v 1.14.6.3 2005/08/17 13:35:52 amoll Exp $
//

#ifndef BALL_STRUCTURE_GRAPHFACE_H
#define BALL_STRUCTURE_GRAPHFACE_H

#ifndef BALL_COMMON_H
#	include <BALL/common.h>
#endif

#include <list>
#include <vector>

namespace BALL
{

	template <typename Vertex, typename Edge, typename Face>
	class GraphVertex;

	template <typename Vertex, typename Edge, typename Face>
	class GraphEdge;

	/** Generic GraphFace Class.
		\ingroup StructureMiscellaneous		
	*/
	template <typename Vertex, typename Edge, typename Face>
	class GraphFace
	{
		public:

		/** @name Class friends

					- class GraphVertex<Vertex,Edge,Face>;
					- class GraphEdge<Vertex,Edge,Face>;
				
		*/
		friend class GraphVertex<Vertex,Edge,Face>;
		friend class GraphEdge<Vertex,Edge,Face>;

		BALL_CREATE(GraphFace)

		/**	@name	Type definitions
		*/
		//@{

		typedef typename std::list<Vertex*>::iterator VertexIterator;
		typedef typename std::list<Vertex*>::const_iterator ConstVertexIterator;
		typedef typename std::list<Edge*>::iterator EdgeIterator;
		typedef typename std::list<Edge*>::const_iterator ConstEdgeIterator;

		//@}
		/**	@name	Constructors and Destructors
		*/
		//@{

		/**	Default constructor.
				This method creates a new GraphFace object.
		*/
		GraphFace()
			throw();

		/**	Copy constructor.
				Create a new GraphFace object from another.
				@param	face	the GraphFace object to be copied
				@param	deep	if deep = false, all pointers are set to NULL (default).
											Otherwise the new GraphFace object is linked to the
											neighbours of the old GraphFace object.
		*/
		GraphFace(const GraphFace<Vertex,Edge,Face>& face, bool deep = false)
			throw();

		/**	Destructor.
				Destructs the GraphFace object.
		*/
		virtual ~GraphFace()
			throw();

		//@}
		/**	@name	Assignments
		*/
		//@{

		/**	Assign from another GraphFace.
				@param	face	the GraphFacee object to assign from
				@param	deep	if deep = false, all pointers are set to NULL	
											(default). Otherwise the GraphFace object is linked to	
											the neighbours of the GraphFace object to assign from.
		*/
		void set(const GraphFace<Vertex,Edge,Face>& face, bool deep = false)
			throw();

		/**	Assign from another GraphFace.
				The GraphFace object is linked to the neighbours of the GraphFace	
				object to assign from.
				@param	face	the GraphFace object to assign from
		*/
		GraphFace<Vertex,Edge,Face>& operator =
				(const GraphFace<Vertex,Edge,Face>& face)
			throw();

		//@}
		/**	@name	Accessors
		*/
		//@{

		/** Insert a new vertex to the GraphFace.
				@param	vertex	a pointer to the new vertex
		*/
		void insert(Vertex* vertex)
			throw();

		/** Insert a new edge to the GraphFace.
				@param	edge	a pointer to the new edge
		*/
		void insert(Edge* edge)
			throw();

		/** Remove a vertex from the GraphFace.
				@param	vertex	a pointer to the vertex to remove
		*/
		void remove(Vertex* vertex)
			throw();

		/** Remove an edge from the GraphFace.
				@param	edge	a pointer to the edge to remove
		*/
		void remove(Edge* edge)
			throw();

		/** Return the number of vertices of the GraphFace.
				@return	Position	the number of vertices of the GraphFace
		*/
		Position numberOfVertices() const
			throw();

		/** Return the number of edges of the GraphFace.
				@return	Position	the number of edges of the GraphFace
		*/
		Position numberOfEdges() const
			throw();


		/** Set the index of the GraphFace.
				@param	index	the new index
		*/
		void setIndex(Index index)
			throw();

		/** Return the index of the GraphFace.
				@return	Index	the index of the GraphFace
		*/
		Index getIndex() const
			throw();

		/** Find the two edges of the GraphFace that belong to the given vertex
				@param	vertex	a pointer to the given vertex
				@param	edge1		a pointer to the first found edge
				@param	edge2		a pointer to the second found edge
				@return	bool		<b>true</b> if the edges can be found,
												<b>false</b> otherwise
		*/
		bool getEdges(const Vertex* vertex, Edge*& edge1, Edge*& edge2) const
			throw();

		/** Find the edge of the GraphFace that belongs to the two given vertices
				@param	vertex1	a pointer to the first given vertex
				@param	vertex2	a pointer to the second given vertex
				@param	edge		a pointer to the found edge
				@return	bool		<b>true</b> if the edge can be found,
												<b>false</b> otherwise
		*/
		bool getEdge
				(const Vertex* vertex1,
				 const Vertex* vertex2,
				 Edge*& edge) const
			throw();

		/** Find the edge of the GraphFace that is similar to the given edge.
				@param	edge	a pointer to the given edge
				@return	Edge*	a pointer to the similar edge of te GraphFace if it can
											be found, otherwise NULL
		*/
		Edge* getSimilarEdge(const Edge* edge) const
			throw();

		/** Substitute a vertex by an other one.
				@param	old_vertex	the vertex that has to be substituted
				@param	new_vertex	the new vertex
				@return	bool				<b>true</b>, if the vertex can be substituted,
														<b>false</b> otherwise
		*/
		bool substitute(const Vertex* old_vertex, Vertex* new_vertex)
			throw();

		/** Substitute an edge by an other one.
				@param	old_edge	the edge that has to be substituted
				@param	new_edge	the new edge
				@return	bool			<b>true</b>, if the edge can be substituted,
													<b>false</b> otherwise
		*/
		bool substitute(const Edge* old_edge, Edge* new_edge)
			throw();

		//@}
		/**	@name	Predicates
		*/
		//@{

		/**	Equality operator.
				@return bool	<b>true</b>
		*/
		virtual bool operator == (const Face& face) const
			throw();

		/**	Inequality operator.
				@return bool	<b>false</b>
		*/
		virtual bool operator != (const Face& face) const
			throw();

		/**	Similarity operator.
				@return bool	<b>true</b>
		*/
		virtual bool operator *= (const Face& face) const
			throw();

		/**	Test whether a vertex is meber of the face.
				@param	vertex	a pointer to the the vertex to test
				@return	Vertex*	a pointer to the vertex if it exists,	
												otherwise <b>NULL</b>
		*/
		Vertex* has(Vertex* vertex) const
			throw();

		/**	Test whether an edge is meber of the face.
				@param	edge	a pointer to the edge to test
				@return	Edge*	a pointer to the edge if it exists, otherwise <b>NULL</b>
		*/
		Edge* has(Edge* edge) const
			throw();

		//@}
		/**	@name	External Iterators
		*/
		//@{

		VertexIterator beginVertex()
			throw();
		ConstVertexIterator beginVertex() const
			throw();
		VertexIterator endVertex()
			throw();
		ConstVertexIterator endVertex() const
			throw();
		EdgeIterator beginEdge()
			throw();
		ConstEdgeIterator beginEdge() const
			throw();
		EdgeIterator endEdge()
			throw();
		ConstEdgeIterator endEdge() const
			throw();

		//@}

		protected:

		/*_	@name	Attributes
		*/
		//@{

		/*_ The vertices of the GraphFace
		*/
		std::list<Vertex*> vertex_;
		/*_ The edges of the GraphFace
		*/
		std::list<Edge*> edge_;
		/* The index of the GraphFace
		*/
		Index index_;

		//@}

	};



	template <typename Vertex, typename Edge, typename Face>
	GraphFace<Vertex,Edge,Face>::GraphFace()
		throw()
		: vertex_(),
			edge_(),
			index_(-1)
	{
	}


	template <typename Vertex, typename Edge, typename Face>
	GraphFace<Vertex,Edge,Face>::GraphFace
			(const GraphFace<Vertex,Edge,Face>& face, bool deep)
		throw()
		: vertex_(),
			edge_(),
			index_(face.index_)
	{
		if (deep)
		{
			vertex_ = face.vertex_;
			edge_ = face.edge_;
		}
	}


	template <typename Vertex, typename Edge, typename Face>
	GraphFace<Vertex,Edge,Face>::~GraphFace()
		throw()
	{
	}


	template <typename Vertex, typename Edge, typename Face>
	void GraphFace<Vertex,Edge,Face>::set
			(const GraphFace<Vertex,Edge,Face>& face, bool deep)
		throw()
	{
		if (this != &face)
		{
			if (deep)
			{
				vertex_ = face.vertex_;
				edge_ = face.edge_;
			}
			index_ = face.index_;
		}
	}


	template <typename Vertex, typename Edge, typename Face>
	GraphFace<Vertex,Edge,Face>& GraphFace<Vertex,Edge,Face>::operator =
			(const GraphFace<Vertex,Edge,Face>& face)
		throw()
	{
		if (this != &face)
		{
			vertex_ = face.vertex_;
			edge_ = face.edge_;
			index_ = face.index_;
		}
		return *this;
	}


	template <typename Vertex, typename Edge, typename Face>
	void GraphFace<Vertex,Edge,Face>::insert(Vertex* vertex)
		throw()
	{
		typename std::list<Vertex*>::iterator v = vertex_.begin();
		bool found = false;
		while ((!found) && (v != vertex_.end()))
		{
			found = (*v == vertex);
			v++;
		}
		if (!found)
		{
			vertex_.push_back(vertex);
		}
	}


	template <typename Vertex, typename Edge, typename Face>
	void GraphFace<Vertex,Edge,Face>::insert(Edge* edge)
		throw()
	{
		typename std::list<Edge*>::iterator e = edge_.begin();
		bool found = false;
		while ((!found) && (e != edge_.end()))
		{
			found = (*e == edge);
			e++;
		}
		if (!found)
		{
			edge_.push_back(edge);
		}
	}


	template <typename Vertex, typename Edge, typename Face>
	void GraphFace<Vertex,Edge,Face>::remove(Vertex* vertex)
		throw()
	{
		vertex_.remove(vertex);
	}


	template <typename Vertex, typename Edge, typename Face>
	void GraphFace<Vertex,Edge,Face>::remove(Edge* edge)
		throw()
	{
		edge_.remove(edge);
	}


	template <typename Vertex, typename Edge, typename Face>
	Position GraphFace<Vertex,Edge,Face>::numberOfVertices() const
		throw()
	{
		return vertex_.size();
	}


	template <typename Vertex, typename Edge, typename Face>
	Position GraphFace<Vertex,Edge,Face>::numberOfEdges() const
		throw()
	{
		return edge_.size();
	}


	template <typename Vertex, typename Edge, typename Face>
	void GraphFace<Vertex,Edge,Face>::setIndex(Index index)
		throw()
	{
		index_ = index;
	}


	template <typename Vertex, typename Edge, typename Face>
	Index GraphFace<Vertex,Edge,Face>::getIndex() const
		throw()
	{
		return index_;
	}


	template <typename Vertex, typename Edge, typename Face>
	bool GraphFace<Vertex,Edge,Face>::getEdges
			(const Vertex*	vertex,
			 Edge*&					edge1,
			 Edge*&					edge2) const
		throw()
	{
		bool found1 = false;
		bool found2 = false;
		typename std::list<Edge*>::const_iterator e = edge_.begin();
		while ((!found1) && (e != edge_.end()))
		{
			if (((*e)->vertex_[0] == vertex) || ((*e)->vertex_[1] == vertex))
			{
				edge1 = *e;
				found1 = true;
			}
			e++;
		}
		if (found1)
		{
			while ((!found2) && (e != edge_.end()))
			{
				if (((*e)->vertex_[0] == vertex) || ((*e)->vertex_[1] == vertex))
				{
					edge2 = *e;
					found2 = true;
				}
				e++;
			}
		}
		return (found1 && found2);
	}


	template <typename Vertex, typename Edge, typename Face>
	bool GraphFace<Vertex,Edge,Face>::getEdge
			(const Vertex*	vertex1,
			 const Vertex*	vertex2,
			 Edge*&					edge) const
		throw()
	{
		typename std::list<Edge*>::const_iterator e = edge_.begin();
		bool found = false;
		while ((!found) && (e != edge_.end()))
		{
			if ((((*e)->vertex_[0] == vertex1) && ((*e)->vertex_[1] == vertex2)) ||
					(((*e)->vertex_[0] == vertex2) && ((*e)->vertex_[1] == vertex1))		)
			{
				edge = *e;
				found = true;
			}
			e++;
		}
		return found;
	}


	template <typename Vertex, typename Edge, typename Face>
	Edge* GraphFace<Vertex,Edge,Face>::getSimilarEdge(const Edge* edge) const
		throw()
	{
		typename std::list<Edge*>::const_iterator e = edge_.begin();
		while (e != edge_.end())
		{
			if (**e *= *edge)
			{
				return *e;
			}
			e++;
		}
		return NULL;
	}


	template <typename Vertex, typename Edge, typename Face>
	bool GraphFace<Vertex,Edge,Face>::substitute
			(const Vertex* old_vertex, Vertex* new_vertex)
		throw()
	{
		typename std::list<Vertex*>::iterator v = vertex_.begin();
		while (v != vertex_.end())
		{
			if (*v == old_vertex)
			{
				*v = new_vertex;
				return true;
			}
			v++;
		}
		return false;
	}


	template <typename Vertex, typename Edge, typename Face>
	bool GraphFace<Vertex,Edge,Face>::substitute
			(const Edge* old_edge, Edge* new_edge)
		throw()
	{
		typename std::list<Edge*>::iterator e = edge_.begin();
		while (e != edge_.end())
		{
			if (*e == old_edge)
			{
				*e = new_edge;
				return true;
			}
			e++;
		}
		return false;
	}


	template <typename Vertex, typename Edge, typename Face>
	bool GraphFace<Vertex,Edge,Face>::operator == (const Face&) const
		throw()
	{
		return true;
	}


	template <typename Vertex, typename Edge, typename Face>
	bool GraphFace<Vertex,Edge,Face>::operator != (const Face&) const
		throw()
	{
		return false;
	}


	template <typename Vertex, typename Edge, typename Face>
	bool GraphFace<Vertex,Edge,Face>::operator *= (const Face&) const
		throw()
	{
		return true;
	}


	template <typename Vertex, typename Edge, typename Face>
	Vertex* GraphFace<Vertex,Edge,Face>::has(Vertex* vertex) const
		throw()
	{
		typename std::list<Vertex*>::const_iterator v = vertex_.begin();
		while (v != vertex_.end())
		{
			if (*v == vertex)
			{
				return *v;
			}
			v++;
		}
		return NULL;
	}


	template <typename Vertex, typename Edge, typename Face>
	Edge* GraphFace<Vertex,Edge,Face>::has(Edge* edge) const
		throw()
	{
		typename std::list<Edge*>::const_iterator e = edge_.begin();
		while (e != edge_.end())
		{
			if (*e == edge)
			{
				return *e;
			}
			e++;
		}
		return NULL;
	}


	template <typename Vertex, typename Edge, typename Face>
	typename GraphFace<Vertex,Edge,Face>::VertexIterator
			GraphFace<Vertex,Edge,Face>::beginVertex()
		throw()
	{
		return vertex_.begin();
	}


	template <typename Vertex, typename Edge, typename Face>
	typename GraphFace<Vertex,Edge,Face>::ConstVertexIterator
			GraphFace<Vertex,Edge,Face>::beginVertex() const
		throw()
	{
		return vertex_.begin();
	}


	template <typename Vertex, typename Edge, typename Face>
	typename GraphFace<Vertex,Edge,Face>::VertexIterator
			GraphFace<Vertex,Edge,Face>::endVertex()
		throw()
	{
		return vertex_.end();
	}


	template <typename Vertex, typename Edge, typename Face>
	typename GraphFace<Vertex,Edge,Face>::ConstVertexIterator
			GraphFace<Vertex,Edge,Face>::endVertex() const
		throw()
	{
		return vertex_.end();
	}


	template <typename Vertex, typename Edge, typename Face>
	typename GraphFace<Vertex,Edge,Face>::EdgeIterator
			GraphFace<Vertex,Edge,Face>::beginEdge()
		throw()
	{
		return edge_.begin();
	}


	template <typename Vertex, typename Edge, typename Face>
	typename GraphFace<Vertex,Edge,Face>::ConstEdgeIterator
			GraphFace<Vertex,Edge,Face>::beginEdge() const
		throw()
	{
		return edge_.begin();
	}


	template <typename Vertex, typename Edge, typename Face>
	typename GraphFace<Vertex,Edge,Face>::EdgeIterator
			GraphFace<Vertex,Edge,Face>::endEdge()
		throw()
	{
		return edge_.end();
	}


	template <typename Vertex, typename Edge, typename Face>
	typename GraphFace<Vertex,Edge,Face>::ConstEdgeIterator
			GraphFace<Vertex,Edge,Face>::endEdge() const
		throw()
	{
		return edge_.end();
	}



	/** Generic GraphTriangle Class.
		\ingroup StructureMiscellaneous		
	*/
	template <typename Vertex, typename Edge, typename Face>
	class GraphTriangle
	{
		public:

		/** @name Class friends

					- class GraphVertex<Vertex,Edge,Face>;
					- class GraphEdge<Vertex,Edge,Face>;
				
		*/
		friend class GraphVertex<Vertex,Edge,Face>;
		friend class GraphEdge<Vertex,Edge,Face>;

		/**	@name	Constructors and Destructors
		*/
		//@{

		BALL_CREATE(GraphTriangle)

		/**	Default constructor.
				This method creates a new GraphTriangle object.
		*/
		GraphTriangle()
			throw();

		/**	Copy constructor.
				Create a new GraphTriangle object from another.
				@param	face	the GraphTriangle object to be copied
				@param	deep	if deep = false, all pointers are set to NULL (default).
											Otherwise the new GraphTriangle object is linked to the
											neighbours of the old GraphTriangle object.
		*/
		GraphTriangle
				(const GraphTriangle<Vertex,Edge,Face>& face, bool deep = false)
			throw();

		/**	Detailled constructor.
				Create a new GraphTriangle object from some nice objects
				@param	vertex1	assigned to the first vertex
				@param	vertex2	assigned to the second vertex
				@param	vertex3	assigned to the third vertex
				@param	edge1		assigned to the first edge
				@param	edge2		assigned to the second edge
				@param	edge3		assigned to the third edge
				@param	index			assigned to the index
		*/
		GraphTriangle
				(Vertex*	vertex1,	Vertex*	vertex2,	Vertex*	vertex3,
				 Edge*		edge1,		Edge*		edge2,		Edge*		edge3,
				 Index		index)
			throw();

		/**	Destructor.
				Destructs the GraphTriangle object.
		*/
		virtual ~GraphTriangle()
			throw();

		//@}
		/**	@name	Assignment
		*/
		//@{

		/**	Assign from another GraphTriangle.
				@param	face	the GraphTriangle object to assign from
				@param	deep	if deep = false, all pointers are set to NULL	
											(default). Otherwise the GraphTriangle object is linked	
											to the neighbours of the GraphTriangle object to assign	
											from.
		*/
		void set(const GraphTriangle<Vertex,Edge,Face>& face, bool deep = false)
			throw();

		/**	Assign from another GraphTriangle.
				The GraphTriangle object is linked to the neighbours of the	
				GraphTriangle object to assifgn from.
				@param	face	the GraphTriangle object to assign from
		*/
		GraphTriangle<Vertex,Edge,Face>& operator =
				(const GraphTriangle<Vertex,Edge,Face>& face)
			throw();

		/**	Assign to a lot of nice objects
				@param	vertex1	assigned to the first vertex
				@param	vertex2	assigned to the second vertex
				@param	vertex3	assigned to the third vertex
				@param	edge1		assigned to the first edge
				@param	edge2		assigned to the second edge
				@param	edge3		assigned to the third edge
				@param	index		assigned to the index
		*/
		void set
				(Vertex*	vertex1,	Vertex*	vertex2,	Vertex*	vertex3,
				 Edge*		edge1,		Edge*		edge2,		Edge*		edge3,
				 Index		index)
			throw();

		//@}
		/**	@name	Accessors
		*/
		//@{

		/** Set one of the vertices of the GraphTriangle.
				@param	i				the relative index of the vertex which should be set.
												If i is greater three, an exception is thrown.
				@param	vertex	a pointer to the new vertex
		*/
		void setVertex(Position i, Vertex* vertex)
			throw(Exception::IndexOverflow);

		/** Return one of the vertices of the GraphTriangle.
				@param	i				the relative index of the vertex which should be given
												back. If i is greater than three, an exception is
												thrown.
				@return	Vertex*	a pointer to the asked vertex
		*/
		Vertex* getVertex(Position i) const
			throw(Exception::IndexOverflow);

		/** Set one of the edges of the GraphTriangle.
				@param	i			the relative index of the edge which should be set.
											If i is greater than three, an exception is thrown.
				@param	edge	a pointer to the new edge
		*/
		void setEdge(Position i, Edge* edge)
			throw(Exception::IndexOverflow);

		/** Return one of the edges of the GraphTriangle.
				@param	i			the relative index of the edge which should be given
											back. If i is greater than three, an exception is thrown.
				@return	Edge*	a pointer to the asked vertex
		*/
		Edge* getEdge(Position i) const
			throw(Exception::IndexOverflow);

		/** Set the index of the GraphTriangle.
				@param	index	the new index
		*/
		void setIndex(Index index)
			throw();

		/** Return the index of the GraphTriangle.
				@return	Index	the index of the GraphTriangle
		*/
		Index getIndex() const
			throw();

		/** Find the two edges of the GraphTriangle that belong to the given vertex.
				@param	vertex	a pointer to the given vertex
				@param	edge1		a pointer to the first found edge
				@param	edge2		a pointer to the second found edge
				@return	bool		<b>true</b> if the edges can be found,
												<b>false</b> otherwise
		*/
		bool getEdges(const Vertex* vertex, Edge*& edge1, Edge*& edge2) const
			throw();

		/** Find the edge of the GraphTriangle that belongs to the two given
				vertices.
				@param	vertex1	a pointer to the first given vertex
				@param	vertex2	a pointer to the second given vertex
				@param	edge		a pointer to the found edge
				@return	bool		<b>true</b> if the edge can be found,
												<b>false</b> otherwise
		*/
		bool getEdge
				(const Vertex* vertex1,
				 const Vertex* vertex2,
				 Edge*& edge) const
			throw();

		/** Find the edge of the GraphFace that is similar to the given edge.
				@param	edge					a pointer to the given edge
				@param	similar_edge	a pointer to the similar edge of te GraphFace
															if it can be found, otherwise NULL
				@return	Index					the relative index of the similar edge if it
															can be found, otherwise -1
		*/
		Index getSimilarEdge(const Edge* edge, Edge*& similar_edge) const
			throw();

		/** Return the relative index of a vertex in the GraphTriangle.
				@return	Index	the relative index of the vertex
		*/
		Index getRelativeIndex(const Vertex* vertex) const
			throw();

		/** Return the relative index of an edge in the GraphTriangle.
				@return	Index	the relative index of the edge
		*/
		Index getRelativeIndex(const Edge* edge) const
			throw();

		/** Return a pointer to the third vertex of the GraphTriangle.
				@param	v1			a pointer to the first vertex
				@param	v2			a pointer to the second vertex
				@return	Vertex*	a pointer to the third vertex
		*/
		Vertex* third(const Vertex* v1, const Vertex* v2) const
			throw();

		/** Return a pointer to the third edge of the GraphTriangle.
				@param	e1		a pointer to the first edge
				@param	e2		a pointer to the second edge
				@return	Edge*	a pointer to the third edge
		*/
		Edge* third(const Edge* e1, const Edge* e2) const
			throw();


		/** Get the edge of the GraphFace which lies on the opposite side of the
				given vertex.
				@param	vertex	a pointer to a vertex of the GraphTriangle
				@return	Edge*		a pointer to the opposite Edge
		*/
		Edge* getOppositeEdge(const Vertex* vertex) const
			throw();


		/** Get the vertex of the GraphFace which lies on the opposite side of the
				given edge.
				@param	edge		a pointer to an edge of the GraphTriangle
				@return	Vertex*	a pointer to the opposite vertex
		*/
		Vertex* getOppositeVertex(const Edge* edge) const
			throw();

		/** Substitute a vertex by an other one.
				@param	old_vertex	the vertex that has to be substituted
				@param	new_vertex	the new vertex
				@return	bool				<b>true</b>, if the vertex can be substituted,
														<b>false</b> otherwise
		*/
		bool substitute(const Vertex* old_vertex, Vertex* new_vertex)
			throw();

		/** Substitute an edge by an other one.
				@param	old_edge	the edge that has to be substituted
				@param	new_edge	the new edge
				@return	bool			<b>true</b>, if the edge can be substituted,
													<b>false</b> otherwise
		*/
		bool substitute(const Edge* old_edge, Edge* new_edge)
			throw();

		//@}
		/**	@name	Predicates
		*/
		//@{

		/**	Equality operator.
				@return bool	<b>true</b>
		*/
		virtual bool operator == (const Face&) const
			throw();

		/**	Inequality operator.
				@return bool	<b>false</b>
		*/
		virtual bool operator != (const Face&) const
			throw();

		/**	Similarity operator.
				@return bool	<b>true</b>
		*/
		virtual bool operator *= (const Face&) const
			throw();

		/**	Test whether a vertex is meber of the face.
				@param	vertex	a pointer to the the vertex to test
				@return	Vertex*	a pointer to the vertex if it exists,	
												otherwise <b>NULL</b>
		*/
		Vertex* has(Vertex* vertex) const
			throw();

		/**	Test whether an edge is meber of the face.
				@param	edge	a pointer to the edge to test
				@return	Edge*	a pointer to the edge if it exists, otherwise <b>NULL</b>
		*/
		Edge* has(Edge* edge) const
			throw();

		//@}

		protected:

		/*_ The vertices of the GraphTriangle
		*/
		Vertex* vertex_[3];
		/*_ The edges of the GraphTriangle
		*/
		Edge* edge_[3];
		/* The index of the GraphTriangle
		*/
		Index index_;

	};



	template <typename Vertex, typename Edge, typename Face>
	GraphTriangle<Vertex,Edge,Face>::GraphTriangle()
		throw()
		: index_(-1)
	{
		vertex_[0] = NULL;
		vertex_[1] = NULL;
		vertex_[2] = NULL;
		edge_[0] = NULL;
		edge_[1] = NULL;
		edge_[2] = NULL;
	}


	template <typename Vertex, typename Edge, typename Face>
	GraphTriangle<Vertex,Edge,Face>::GraphTriangle
			(const GraphTriangle<Vertex,Edge,Face>& face, bool deep)
		throw()
		: index_(face.index_)
	{
		if (deep)
		{
			vertex_[0] = face.vertex_[0];
			vertex_[1] = face.vertex_[1];
			vertex_[2] = face.vertex_[2];
			edge_[0] = face.edge_[0];
			edge_[1] = face.edge_[1];
			edge_[2] = face.edge_[2];
		}
		else
		{
			vertex_[0] = NULL;
			vertex_[1] = NULL;
			vertex_[2] = NULL;
			edge_[0] = NULL;
			edge_[1] = NULL;
			edge_[2] = NULL;
		}
	}


	template <typename Vertex, typename Edge, typename Face>
	GraphTriangle<Vertex,Edge,Face>::GraphTriangle
			(Vertex*	vertex1,	Vertex*	vertex2,	Vertex*	vertex3,
			 Edge*		edge1,		Edge*		edge2,		Edge*		edge3,
			 Index		index)
		throw()
		: index_(index)
	{
		vertex_[0] = vertex1;
		vertex_[1] = vertex2;
		vertex_[2] = vertex3;
		edge_[0] = edge1;
		edge_[1] = edge2;
		edge_[2] = edge3;
	}


	template <typename Vertex, typename Edge, typename Face>
	GraphTriangle<Vertex,Edge,Face>::~GraphTriangle()
		throw()
	{
	}


	template <typename Vertex, typename Edge, typename Face>
	void GraphTriangle<Vertex,Edge,Face>::set
			(const GraphTriangle<Vertex,Edge,Face>& face, bool deep)
		throw()
	{
		if (this != &face)
		{
			if (deep)
			{
				vertex_[0] = face.vertex_[0];
				vertex_[1] = face.vertex_[1];
				vertex_[2] = face.vertex_[2];
				edge_[0] = face.edge_[0];
				edge_[1] = face.edge_[1];
				edge_[2] = face.edge_[2];
			}
			else
			{
				vertex_[0] = NULL;
				vertex_[1] = NULL;
				vertex_[2] = NULL;
				edge_[0] = NULL;
				edge_[1] = NULL;
				edge_[2] = NULL;
			}
			index_ = face.index_;
		}
	}


	template <typename Vertex, typename Edge, typename Face>
	GraphTriangle<Vertex,Edge,Face>& GraphTriangle<Vertex,Edge,Face>::operator =
			(const GraphTriangle<Vertex,Edge,Face>& face)
		throw()
	{
		if (this != &face)
		{
			vertex_[0] = face.vertex_[0];
			vertex_[1] = face.vertex_[1];
			vertex_[2] = face.vertex_[2];
			edge_[0] = face.edge_[0];
			edge_[1] = face.edge_[1];
			edge_[2] = face.edge_[2];
			index_ = face.index_;
		}
		return *this;
	}


	template <typename Vertex, typename Edge, typename Face>
	void GraphTriangle<Vertex,Edge,Face>::set
			(Vertex*	vertex1,	Vertex*	vertex2, Vertex*	vertex3,
			 Edge*		edge1,		Edge*		edge2,	 Edge*		edge3,
			 Index		index)
		throw()
	{
		vertex_[0] = vertex1;
		vertex_[1] = vertex2;
		vertex_[2] = vertex3;
		edge_[0] = edge1;
		edge_[1] = edge2;
		edge_[2] = edge3;
		index_ = index;
	}


	template <typename Vertex, typename Edge, typename Face>
	void GraphTriangle<Vertex,Edge,Face>::setVertex(Position i, Vertex* vertex)
		throw(Exception::IndexOverflow)
	{
		if (i > 2)
		{
			throw Exception::IndexOverflow(__FILE__,__LINE__,i,2);
		}
		else
		{
			vertex_[i] = vertex;
		}
	}


	template <typename Vertex, typename Edge, typename Face>
	Vertex* GraphTriangle<Vertex,Edge,Face>::getVertex(Position i) const
		throw(Exception::IndexOverflow)
	{
		if (i > 2)
		{
			throw Exception::IndexOverflow(__FILE__,__LINE__,i,2);
		}
		else
		{
			return vertex_[i];
		}
	}


	template <typename Vertex, typename Edge, typename Face>
	void GraphTriangle<Vertex,Edge,Face>::setEdge(Position i, Edge* edge)
		throw(Exception::IndexOverflow)
	{
		if (i > 2)
		{
			throw Exception::IndexOverflow(__FILE__,__LINE__,i,2);
		}
		else
		{
			edge_[i] = edge;
		}
	}


	template <typename Vertex, typename Edge, typename Face>
	Edge* GraphTriangle<Vertex,Edge,Face>::getEdge(Position i) const
		throw(Exception::IndexOverflow)
	{
		if (i > 2)
		{
			throw Exception::IndexOverflow(__FILE__,__LINE__,i,2);
		}
		else
		{
			return edge_[i];
		}
	}


	template <typename Vertex, typename Edge, typename Face>
	void GraphTriangle<Vertex,Edge,Face>::setIndex(Index index)
		throw()
	{
		index_ = index;
	}


	template <typename Vertex, typename Edge, typename Face>
	Index GraphTriangle<Vertex,Edge,Face>::getIndex() const
		throw()
	{
		return index_;
	}


	template <typename Vertex, typename Edge, typename Face>
	bool GraphTriangle<Vertex,Edge,Face>::getEdges
			(const Vertex*	vertex,
			 Edge*&					edge1,
			 Edge*&					edge2) const
		throw()
	{
		Position i = 0;
		bool found1 = false;
		bool found2 = false;
		while ((!found1) && (i < 3))
		{
			if (edge_[i] != NULL)
			{
				if ((edge_[i]->vertex_[0] == vertex) ||
						(edge_[i]->vertex_[1] == vertex)		)
				{
					edge1 = edge_[i];
					found1 = true;
				}
			}
			i++;
		}
		if (found1)
		{
			while ((!found2) && (i < 3))
			{
				if (edge_[i] != NULL)
				{
					if ((edge_[i]->vertex_[0] == vertex) ||
							(edge_[i]->vertex_[1] == vertex)		)
					{
						edge2 = edge_[i];
						found2 = true;
					}
				}
				i++;
			}
		}
		return (found1 && found2);
	}


	template <typename Vertex, typename Edge, typename Face>
	bool GraphTriangle<Vertex,Edge,Face>::getEdge
			(const Vertex*	vertex1,
			 const Vertex*	vertex2,
			 Edge*&					edge) const
		throw()
	{
		Position i = 0;
		bool found = false;
		while ((!found) && (i < 3))
		{
			if (edge_[i] != NULL)
			{
				if (((edge_[i]->vertex_[0] == vertex1) &&
						 (edge_[i]->vertex_[1] == vertex2)		) ||
						((edge_[i]->vertex_[0] == vertex2) &&
						 (edge_[i]->vertex_[1] == vertex1)		)			)
				{
					edge = edge_[i];
					found = true;
				}
			}
			i++;
		}
		return found;
	}


	template <typename Vertex, typename Edge, typename Face>
	Index GraphTriangle<Vertex,Edge,Face>::getSimilarEdge
			(const Edge* edge, Edge*& similar_edge) const
		throw()
	{
		if (*edge_[0] *= *edge)
		{
			similar_edge = edge_[0];
			return 0;
		}
		if (*edge_[1] *= *edge)
		{
			similar_edge = edge_[1];
			return 1;
		}
		if (*edge_[2] *= *edge)
		{
			similar_edge = edge_[2];
			return 2;
		}
		similar_edge = NULL;
		return -1;
	}


	template <typename Vertex, typename Edge, typename Face>
	Index GraphTriangle<Vertex,Edge,Face>::getRelativeIndex
			(const Vertex* vertex) const
		throw()
	{
		for (Position i = 0; i < 3; i++)
		{
			if (vertex_[i] == vertex)
			{
				return i;
			}
		}
		return -1;
	}


	template <typename Vertex, typename Edge, typename Face>
	Index GraphTriangle<Vertex,Edge,Face>::getRelativeIndex
			(const Edge* edge) const
		throw()
	{
		for (Position i = 0; i < 3; i++)
		{
			if (edge_[i] == edge)
			{
				return i;
			}
		}
		return -1;
	}


	template <typename Vertex, typename Edge, typename Face>
	Vertex* GraphTriangle<Vertex,Edge,Face>::third
			(const Vertex* v1, const Vertex* v2) const
		throw()
	{
		if ((vertex_[0] == v1) || (vertex_[0] == v2))
		{
			if ((vertex_[1] == v1) || (vertex_[1] == v2))
			{
				return vertex_[2];
			}
			else
			{
				return vertex_[1];
			}
		}
		else
		{
			return vertex_[0];
		}
	}


	template <typename Vertex, typename Edge, typename Face>
	Edge* GraphTriangle<Vertex,Edge,Face>::third
			(const Edge* e1, const Edge* e2) const
		throw()
	{
		if ((edge_[0] == e1) || (edge_[0] == e2))
		{
			if ((edge_[1] == e1) || (edge_[1] == e2))
			{
				return edge_[2];
			}
			else
			{
				return edge_[1];
			}
		}
		else
		{
			return edge_[0];
		}
	}


	template <typename Vertex, typename Edge, typename Face>
	Edge* GraphTriangle<Vertex,Edge,Face>::getOppositeEdge
			(const Vertex* vertex) const
		throw()
	{
		for (Position i = 0; i < 3; i++)
		{
			if ((edge_[i]->vertex_[0] != vertex) &&
					(edge_[i]->vertex_[1] != vertex)		)
			{
				return edge_[i];
			}
		}
		return NULL;
	}


	template <typename Vertex, typename Edge, typename Face>
	Vertex* GraphTriangle<Vertex,Edge,Face>::getOppositeVertex
			(const Edge* edge) const
		throw()
	{
		for (Position i = 0; i < 3; i++)
		{
			if ((vertex_[i] != edge->vertex_[0]) &&
					(vertex_[i] != edge->vertex_[1])		)
			{
				return vertex_[i];
			}
		}
		return NULL;
	}


	template <typename Vertex, typename Edge, typename Face>
	bool GraphTriangle<Vertex,Edge,Face>::substitute
			(const Vertex* old_vertex, Vertex* new_vertex)
		throw()
	{
		for (Position i = 0; i < 3; i++)
		{
			if (vertex_[i] == old_vertex)
			{
				vertex_[i] = new_vertex;
				return true;
			}
		}
		return false;
	}


	template <typename Vertex, typename Edge, typename Face>
	bool GraphTriangle<Vertex,Edge,Face>::substitute
			(const Edge* old_edge, Edge* new_edge)
		throw()
	{
		for (Position i = 0; i < 3; i++)
		{
			if (edge_[i] == old_edge)
			{
				edge_[i] = new_edge;
				return true;
			}
		}
		return false;
	}


	template <typename Vertex, typename Edge, typename Face>
	bool GraphTriangle<Vertex,Edge,Face>::operator == (const Face&) const
		throw()
	{
		return true;
	}


	template <typename Vertex, typename Edge, typename Face>
	bool GraphTriangle<Vertex,Edge,Face>::operator != (const Face&) const
		throw()
	{
		return false;
	}


	template <typename Vertex, typename Edge, typename Face>
	bool GraphTriangle<Vertex,Edge,Face>::operator *= (const Face&) const
		throw()
	{
		return true;
	}


	template <typename Vertex, typename Edge, typename Face>
	Vertex* GraphTriangle<Vertex,Edge,Face>::has(Vertex* vertex) const
		throw()
	{
		if (vertex_[0] == vertex)
		{
			return vertex_[0];
		}
		if (vertex_[1] == vertex)
		{
			return vertex_[1];
		}
		if (vertex_[2] == vertex)
		{
			return vertex_[2];
		}
		return NULL;
	}


	template <typename Vertex, typename Edge, typename Face>
	Edge* GraphTriangle<Vertex,Edge,Face>::has(Edge* edge) const
		throw()
	{
		if (edge_[0] == edge)
		{
			return edge_[0];
		}
		if (edge_[1] == edge)
		{
			return edge_[1];
		}
		if (edge_[2] == edge)
		{
			return edge_[2];
		}
		return NULL;
	}


} // namespace BALL

#endif // BALL_STRUCTURE_RSFACE_H
