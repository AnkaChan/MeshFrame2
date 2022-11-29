/*!
*      \file halfface.h
*      \brief Base HalfFace Class for all types of Tetrahedron Mesh Classes
*
*		This is the fundamental class for HalfFace
*	   \author David Gu
*      \date 10/01/2011
*
*/

#ifndef _TMESHLIB_HALFFACE_H_
#define _TMESHLIB_HALFFACE_H_

#include <list>

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
		* \brief CHalfFace, base class for HalfFace
		*/
		class CHalfFace
		{
		public:
			CHalfFace()
			{
				m_pHalfEdge = nullptr;
				m_pFace = nullptr;
				m_pTet = nullptr;
				m_pDual = nullptr;
			};

			~CHalfFace(){};

			CHalfEdge * half_edge() { return m_pHalfEdge; };
			CFace     * face()      { return m_pFace; };
			CTet      * tet()       { return m_pTet; };
			CHalfFace * dual()      { return m_pDual; };
			int       & key(int k) { return m_key[k]; };

			void SetHalfEdge(CHalfEdge * pHe) { m_pHalfEdge = pHe; };
			void SetFace(CFace * pF)          { m_pFace = pF; };
			void SetTet(CTet  * pT)          { m_pTet = pT; };
			void SetDual(CHalfFace   *pF)     { m_pDual = pF; };
			bool operator==(const CHalfFace & f)
			{
				for (int i = 0; i < 3; i++)
					if (m_key[i] != f.m_key[i]) return false;
				return true;
			};

			void setSurfaceMeshFacePtr(void * in_pSurfaceMeshFace) {
				m_surfaceFace = in_pSurfaceMeshFace;
			}

			void* getSurfaceFacePtr() {
				return m_surfaceFace;
			}

		public:
			size_t& index() { return m_index; };
		private:
			size_t m_index;


		protected:

			CHalfEdge * m_pHalfEdge;
			CFace     * m_pFace;
			CTet      * m_pTet;
			CHalfFace * m_pDual;
			int         m_key[3];
			void      * m_surfaceFace = nullptr;
		};
	};
};

#endif