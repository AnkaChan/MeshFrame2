#include "BaseTMesh.h"
#include "Eigen/core"

#ifndef _MESHFRAME_STATIC_TET_MESH_H_
#define _MESHFRAME_STATIC_TET_MESH_H_

namespace MF
{
	namespace TetMesh
	{

		template <typename DType, typename TVertexType, typename VertexType, typename HalfEdgeType, typename TEdgeType, typename EdgeType, typename HalfFaceType, typename FaceType, typename TetType>
		class CTMeshStatic : public class CTMeshBase<TVertexType, VertexType, HalfEdgeType, TEdgeType, EdgeType, HalfFaceType, FaceType, TetType>
		{
		public:
			typedef CTMeshStatic<TVType, VType, HEType, TEType, EType, HFType, FType, TType>* Ptr;
			typedef std::shared_ptr<CTMeshStatic<TVType, VType, HEType, TEType, EType, HFType, FType, TType>> SharedPtr;


			typedef Eigen::Matrix<DType, POINT_VEC_DIMS, Eigen::Dynamic> TVerticesMat;


		private:
			TVerticesMat mVertPos;


		};
	}
}

#endif