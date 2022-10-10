/*!
*      \file tvertex.h
*      \brief Base tvertex Class for all types of Tetrahedron Mesh Classes
*
*		This is the fundamental class for tvertex
*	   \author David Gu
*      \date 10/01/2011
*
*/

#ifndef _TMESHLIB_TVERTEX_H_
#define _TMESHLIB_TVERTEX_H_

#include <list>

namespace MF
{
	namespace TetMesh
	{

		class CVertexBase;
		class CHalfEdge;
		class CEdge;
		class CHalfFace;
		class CFace;
		class CTet;


		/*!
		* \brief CTVertex, base class for Tetrahedron vertex
		*/
		class CTVertex
		{
		public:
			CTVertex() { m_pVertex = NULL; m_pTet = NULL; m_pHalfedge = NULL; };
			~CTVertex() {};

			int & id() { return m_iID; };

			CVertexBase* vert() { return m_pVertex; };
			CTet            *  tet() { return m_pTet; };
			CHalfEdge       * halfedge() { return m_pHalfedge; };

			void set_vert(CVertexBase * pV) { m_pVertex = pV; };
			void set_tet(CTet        * pT) { m_pTet = pT; };
			void set_halfedge(CHalfEdge * pH) { m_pHalfedge = pH; };

		public:
			size_t& index() { return m_index; };
		private:
			size_t m_index;
		protected:
			//vertex ID
			int            m_iID;
			CVertexBase* m_pVertex;
			CTet         * m_pTet;
			//outgoing, halfedge start from this TVertex
			CHalfEdge *    m_pHalfedge;
		};

	};
};

#endif