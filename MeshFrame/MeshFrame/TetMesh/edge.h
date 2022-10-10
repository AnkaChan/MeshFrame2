/*!
*      \file edge.h
*      \brief Base edge Class for all types of Tetrahedron Mesh Classes
*
*		This is the fundamental class for edge
*	   \author David Gu
*      \date 10/01/2011
*
*/

#ifndef _TMESHLIB_EDGE_H_
#define _TMESHLIB_EDGE_H_

#include <list>
#include <vector>


namespace MF
{
	namespace TetMesh
	{
		class CVertexBase;
		class CTVertex;
		class CHalfEdge;
		class CEdge;
		class CHalfFace;
		class CFace;
		class CTet;

		/*!
		* \brief CEdge, base class for edge
		*/
		class CEdge	//Edge among tets
		{
		public:

			CEdge() { m_vertices[0] = NULL; m_vertices[1] = NULL;  m_bIsBoundary = false; };

			~CEdge() { m_lTEdges.clear(); };

			std::vector<CTEdge*> * edges() { return &m_lTEdges; };

			bool   & boundary() { return m_bIsBoundary; };

			CVertexBase*  vertex1() { return m_vertices[0]; };
			CVertexBase*  vertex2() { return m_vertices[1]; };

			void SetVertex1(CVertexBase* v) { m_vertices[0] = v; };
			void SetVertex2(CVertexBase* v) { m_vertices[1] = v; };

		public:
			size_t& index() { return m_index; };
		private:
			size_t m_index;
		protected:
			std::vector<CTEdge*> m_lTEdges;
			CVertexBase*  m_vertices[2];
			bool     m_bIsBoundary;

		};

	};
};

#endif