#ifndef _TMESHLIB_HEADERS_H_
#define _TMESHLIB_HEADERS_H_

#include "TMeshStatic.h"
#include "VertexStatic.h"
#include "Tvertex.h"
#include "Edge.h"
#include "TEdge.h"
#include "Face.h"
#include "Halfface.h"
#include "Halfedge.h"
#include "Tet.h"

#include "TIterators.h"

namespace MF
{
	namespace TetMesh
	{
		template <typename DType>
		using CTMeshStaticDType =  CTMeshStatic<DType, CTVertex, CVertexStatic<DType>, CHalfEdge, CTEdge, CEdge, CHalfFace, CFace, CTet>;
		typedef CTMeshStaticDType<float> TMeshStaticF;
		typedef CTMeshStaticDType<double> TMeshStaticD;
	}
}

#endif