#include <Eigen/core>
#include <Eigen/Dense>

#define POINT_VEC_DIMS 3


namespace MF
{
	template<typename DType>
	using TVerticesMat = Eigen::Matrix<DType, POINT_VEC_DIMS, Eigen::Dynamic>;

	template<typename DType>
	using TVec3 = Eigen::Matrix<DType, 3, 1>;
	template<typename DType>
	using TVec2 = Eigen::Matrix<DType, 2, 1>;

	template<typename DType>
	using TVec3Block = Eigen::Block<TVerticesMat<DType>, POINT_VEC_DIMS, 1>;

}