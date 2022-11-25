#ifndef _MESHFRAME_MESHSTATIC_H_
#define _MESHFRAME_MESHSTATIC_H_

#include "BaseMesh.h"

namespace MF {
	namespace TriMesh {
		// VertexType mush be child class of VertexStatic
		template <typename DType, typename VertexType, typename EdgeType, typename FaceType, typename HalfEdgeType>
		class CMeshStatic : CMeshBase<DType, VertexType, EdgeType, FaceType, HalfEdgeType> {

		};
	}

}

#endif