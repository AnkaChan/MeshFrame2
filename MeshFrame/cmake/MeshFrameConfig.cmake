message( "Adding MeshFrame." )

# please use x64 arch
SET (MESHFRAME_INCLUDE_DIR  ${CMAKE_CURRENT_LIST_DIR}/..)


file(GLOB MESHFRAME_SOURCE_CPP_UTILITY
    ${MESHFRAME_INCLUDE_DIR}/MeshFrame/Utility/*.cpp
)

SET (MESHFRAME_SOURCE_CPP
	${MESHFRAME_SOURCE_CPP_UTILITY}
)

# SET (MESHFRAME_SOURCE_CPP_CUDA
# 	${MESHFRAME_SOURCE_CPP}/MeshFrame/common/logger.cpp
# 	${MESHFRAME_SOURCE_CPP}/MeshFrame/common/windows/getopt.c
# 	${MESHFRAME_SOURCE_CPP}/MeshFrame/NetWrapper/NetWrapper.cpp
# 	${MESHFRAME_SOURCE_CPP}/MeshFrame/MarkerLabeler/MarkerLabeler.cpp
# 	${MESHFRAME_SOURCE_CPP}/MeshFrame/Image/Image.cpp
# 	
# )

find_package(Eigen3 REQUIRED)

message( "Found Eigen3 at:"  ${EIGEN3_INCLUDE_DIRS})

SET (MESHFRAME_INCLUDE_DIR 
	${MESHFRAME_INCLUDE_DIR}
	${EIGEN3_INCLUDE_DIRS}
)





