/*!
*      \file vertex.h
*      \brief Base vertex Class for all types of Tetrahedron Mesh Classes
*
*		This is the fundamental class for vertex
*	   \date 10/10/2022
*
*/

#ifndef _TMESHLIB_VERTEX_STATIC_H_
#define _TMESHLIB_VERTEX_STATIC_H_

#include <list>
#include <vector>
#include "../Memory/Array.h"
#include "vertex.h"
#include "../Types/TypeDefs.h"

namespace MF
{
	namespace TetMesh
	{
		template <typename DType>
		class CVertexStatic: public CVertexBase
		{
		public:
			CVertexStatic() { };
			~CVertexStatic(){ }

			void setPVertPos(TVerticesMat<DType>* pVertPos) {
				mPVertPos = pVertPos;
			}

			TVec3Block<DType> position() {
				return mPVertPos->block<3, 1>(0, m_id);
			}

			 
		protected:
			
			TVerticesMat<DType>* mPVertPos = nullptr;
		};
	};
};

#endif