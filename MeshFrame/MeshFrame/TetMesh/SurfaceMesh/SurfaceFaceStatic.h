/*!
*      \file Vertex.h
*      \brief Base class of vertex
*	   \Version 2.0
*	   \Update 11/25/2022
*/

#ifndef  _MESHFRAME_TETMESH_SURFACE_FACE_H_
#define  _MESHFRAME_TETMESH_SURFACE_FACE_H_

#include "../../TriMesh/Face.h"
#include "../halfface.h"

namespace MF {
	namespace TetMesh {
		class CSurfaceFace : public TriMesh::CFace
		{
		public:
			CSurfaceFace() { };
			~CSurfaceFace() { }

			void setTetMeshHalfFacePtr(TetMesh::CHalfFace* in_PTetMeshVertex) {
				pTetMeshHalfFace = in_PTetMeshVertex;
				in_PTetMeshVertex->setSurfaceMeshFacePtr((void*)this);
			}
			TetMesh::CHalfFace* getTetMeshHalfFacePtr() {
				return pTetMeshHalfFace;
			}


		protected:
			TetMesh::CHalfFace * pTetMeshHalfFace = nullptr;

		};
	}
}

#endif //_MESHLIB_VERTEX_H_defined