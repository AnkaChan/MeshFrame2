#include "Eigen/core"
#include "../Types/TypeDefs.h"


namespace MF
{
	namespace TetMesh {
		typedef Eigen::Matrix<IdType, 4, Eigen::Dynamic> TTetIdsMat;
		using Vec4BlockI = Eigen::Block<TTetIdsMat, 4, 1>;

	}

}