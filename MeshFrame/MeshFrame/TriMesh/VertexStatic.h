/*!
*      \file Vertex.h
*      \brief Base class of vertex
*	   \Version 2.0
*	   \Update 11/03/2022
*/

#ifndef  _MESHFRAME_VERTEX_STATIC_H_
#define  _MESHFRAME_VERTEX_STATIC_H_

#include "Vertex.h"
#include "../Types/TypeDefs.h"

#include "MeshStatic.h"

namespace MF {
	namespace TriMesh {
		template <typename DType>
		class CVertexStatic : public TriMesh::CVertexBase
		{
		public:
			void setPVertPos(TVerticesMat<DType>* pVertPos) {
				mPVertPos = pVertPos;
			}

			TVec3Block<DType> position() {
				return mPVertPos->block<3, 1>(0, m_id);
			}

			virtual void setData(void* pMesh);

		protected:
			TVerticesMat<DType>* mPVertPos = nullptr;
		};

	}

}//name space MF

#endif //_MESHFRAME_VERTEX_STATIC_H_