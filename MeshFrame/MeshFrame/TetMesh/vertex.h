/*!
*      \file vertex.h
*      \brief Base vertex Class for all types of Tetrahedron Mesh Classes
*
*		This is the fundamental class for vertex
*	   \date 10/10/2022
*
*/

#ifndef _TMESHLIB_VERTEX_H_
#define _TMESHLIB_VERTEX_H_

#include <list>
#include <vector>
#include "../Memory/Array.h"

namespace MF
{
	namespace TetMesh
	{
		class CTVertex;
		class CHalfEdge;
		class CEdge;
		class CTEdge;
		class CHalfFace;
		class CFace;
		class CTet;

		typedef std::vector<CEdge*> CEArray;
		typedef std::vector<CTVertex*> CTVArray;

		class CVertexBase
		{
		public:
			CVertexBase() { m_index = 0; m_bIsBoundary = false; };
			~CVertexBase(){ }

			virtual int& id()       { return m_id; };
			size_t & index()    { return m_index; };
			bool   & boundary()	{ return m_bIsBoundary; };

			CEArray * edges() { return &m_pEdges; };
			CTVArray * tvertices() { return &m_pTVertices; };
			std::list<CHalfFace*> * HalfFaces(){ return &m_pHFaces; };

		protected:
			int    m_id;
			size_t    m_index;
			bool   m_bIsBoundary;

			std::list<CHalfFace*>  m_pHFaces;		//temporary HalfFace list, will be empty after loading the whole mesh 

			CTVArray   m_pTVertices;	//adjacent TVertecies
			CEArray      m_pEdges;	    //adjacent Edges;

			//std::string m_string;
		};
	};
};

#endif