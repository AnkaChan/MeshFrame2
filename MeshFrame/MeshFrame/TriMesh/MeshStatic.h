#ifndef _MESHFRAME_MESHSTATIC_H_
#define _MESHFRAME_MESHSTATIC_H_

#include "BaseMesh.h"
#include "VertexStatic.h"
#include "Face.h"
#include "Edge.h"
#include "HalfEdge.h"

#include "Iterators2.h"

namespace MF {
	namespace TriMesh {
		template <typename DType>
		using CTriMeshStaticDType = CMeshBase<DType, CVertexStatic<DType>, CEdge, CFace, CHalfEdge>;
		typedef CTriMeshStaticDType<float> TriMeshStaticF;
		typedef CTriMeshStaticDType<double> TriMeshStaticD;

		typedef CIterators<TriMeshStaticF>  TriMeshStaticIteratorF;
		typedef CIterators<TriMeshStaticD>  TriMeshStaticIteratorD;
	}
}

#endif