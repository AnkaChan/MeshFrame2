/*!
*      \file Vertex.h
*      \brief Base class of vertex
*	   \Version 2.0
*	   \Update 11/25/2022
*/

#ifndef  _MESHFRAME_TETMESH_SURFACE_VERTEX_H_
#define  _MESHFRAME_TETMESH_SURFACE_VERTEX_H_

#include "../../TriMesh/vertex.h"
#include "../VertexStatic.h"

namespace MF {
	namespace TetMesh {
		template <typename DType>
		class CSurfaceVertexStatic : public TriMesh::CVertexBase
		{
		public:
			CSurfaceVertexStatic() { };
			~CSurfaceVertexStatic() { }

			void setTetMeshVertPtr(TetMesh::CVertexStatic<DType>* in_PTetMeshVertex) {
				pTetMeshVertex = in_PTetMeshVertex;
				in_PTetMeshVertex->setSurfaceMeshVertexPtr((void*)this);
			}

			TVec3Block<DType> position() {
				return pTetMeshVertex->position();
			}

		protected:
			TetMesh::CVertexStatic<DType>* pTetMeshVertex = nullptr;

		};
	}
}

#endif //_MESHLIB_VERTEX_H_defined