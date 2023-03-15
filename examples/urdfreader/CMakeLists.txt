PROJECT (RBDLEXAMPLE CXX)

CMAKE_MINIMUM_REQUIRED(VERSION 3.0)

# We need to add the project source path to the CMake module path so that
# the FindRBDL.cmake script can be found.
LIST( APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR} )


SET(CUSTOM_RBDL_PATH "" CACHE PATH "Path to specific RBDL Installation")
ADD_EXECUTABLE (example_urdfreader example_urdfreader.cc)

# Search for the RBDL include directory and library
FIND_PACKAGE (RBDL COMPONENTS URDFREADER REQUIRED)
FIND_PACKAGE (Eigen3 REQUIRED)
find_package(Python2 COMPONENTS Interpreter Development REQUIRED)
find_package(Python2 COMPONENTS NumPy)
# if(Python3_NumPy_FOUND)
#   target_link_libraries(example_urdfreader
#     Python3::NumPy
#   )
# else()
#   target_compile_definitions(example_urdfreader INTERFACE WITHOUT_NUMPY)
# endif()

# Add the include directory to the include paths
INCLUDE_DIRECTORIES ( ${RBDL_INCLUDE_DIR} ${EIGEN3_INCLUDE_DIR})

# Create an executable
# ADD_EXECUTABLE (example_urdfreader example_urdfreader.cc)

# And link the library against the executable
TARGET_LINK_LIBRARIES ( example_urdfreader
	${RBDL_LIBRARY}
	${RBDL_URDFREADER_LIBRARY}
	Python2::Python
	Python2::Module
	Python2::NumPy
	)
