#ifndef _TMESHLIB_SURFACE_MESH_HEADERS_H_
#define _TMESHLIB_SURFACE_MESH_HEADERS_H_

#include "../TMeshStaticLibHeaders.h"
#include "../../TriMesh/MeshCoreHeaders.h"

#include "SurfaceVertexStatic.h"
#include "SurfaceFaceStatic.h"
#include "SurfaceMeshStatic.h"


namespace MF {
	namespace TetMesh {
		//template <typename DType, typename VertexType, typename EdgeType, typename FaceType, typename HalfEdgeType>
		//class CSurfaceMeshStatic : TriMesh::CMeshBase<DType, VertexType, EdgeType, FaceType, HalfEdgeType>

		template <typename DType>
		using CSurfaceMeshStaticDType = CSurfaceMeshStatic<DType, CSurfaceVertexStatic<DType>, TriMesh::CEdge, CSurfaceFace, TriMesh::CHalfEdge>;

		typedef CSurfaceMeshStaticDType<float> CSurfaceMeshStaticF;
		typedef CSurfaceMeshStaticDType<double> CSurfaceMeshStaticD;
	}
}

#endif