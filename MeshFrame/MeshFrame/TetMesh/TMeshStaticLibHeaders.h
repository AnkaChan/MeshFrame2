#ifndef _TMESHLIB_HEADERS_H_
#define _TMESHLIB_HEADERS_H_

#include "BaseTMesh.h"
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
		typedef CTMeshStatic<float, CTVertex, CVertexStatic<float>, CHalfEdge, CTEdge, CEdge, CHalfFace, CFace, CTet> TMeshStaticF;
		typedef CTMeshStatic<double, CTVertex, CVertexStatic<double>, CHalfEdge, CTEdge, CEdge, CHalfFace, CFace, CTet> TMeshStaticD;
	}
}

#endif