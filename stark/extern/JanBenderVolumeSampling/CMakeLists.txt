################################################################################
# CMake setup
################################################################################
cmake_minimum_required(VERSION 3.9)

project(JanBenderVolumeSampling)

# disable if too many error messages from other libraries
# add_compile_options(-Wall -Wextra -pedantic)

##  Includes and executes CMake/common.cmake
include(CMake/Common.cmake)
if (NOT WIN32)
    message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")
endif()

## Specifies the location of external installations needed by the program
set(ExternalInstallDir "${CMAKE_BINARY_DIR}/extern/install" CACHE INTERNAL "")
set(EXT_CMAKE_BUILD_TYPE ${CMAKE_BUILD_TYPE} CACHE INTERNAL "")
if (NOT ${CMAKE_BUILD_TYPE} STREQUAL "Debug")
    set(EXT_CMAKE_BUILD_TYPE "Release" CACHE INTERNAL "")
endif()

################################################################################
# Foreign external libraries
################################################################################

## Executes FindEigen3.cmake in CMake/ and disables errors due to alignmetn
find_package(Eigen3 REQUIRED)
add_definitions(-DEIGEN_DISABLE_UNALIGNED_ARRAY_ASSERT)

################################################################################
# Own external libraries (Executes CMake/SetUpExternalProjects.cmake)
################################################################################
include(CMake/SetUpExternalProjects.cmake)

set(JBVS_LINK_LIBRARIES ${DISCREGRID_LIBRARIES})
set(JBVS_DEPENDENCIES Ext_Discregrid)

set(JBVS_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/include)
set(JBVS_SRC_DIR ${CMAKE_CURRENT_SOURCE_DIR}/src)

set(JBVS_INCLUDE_FILES
    ${JBVS_INCLUDE_DIR}/SDFFunctions.h
    ${JBVS_INCLUDE_DIR}/VolumeSampling.h
)

set(JBVS_SRC_FILES
    ${JBVS_SRC_DIR}/SDFFunctions.cpp
    ${JBVS_SRC_DIR}/VolumeSampling.cpp
)

add_library(JanBenderVolumeSampling
    ${JBVS_INCLUDE_FILES}
    ${JBVS_SRC_FILES}
)

add_dependencies(JanBenderVolumeSampling ${JBVS_DEPENDENCIES})

target_link_libraries(JanBenderVolumeSampling PUBLIC ${JBVS_LINK_LIBRARIES})

target_include_directories(JanBenderVolumeSampling PUBLIC
        ${EIGEN3_INCLUDE_DIR}
        ${DISCREGRID_INCLUDE_DIR}
        ${CMAKE_SOURCE_DIR}/stark/extern/JanBenderUtilities/include
        ${JBVS_INCLUDE_DIR}
)