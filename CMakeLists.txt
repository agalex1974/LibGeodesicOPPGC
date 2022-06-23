cmake_minimum_required(VERSION 3.22)

set(default_build_type "Release")

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    message(STATUS "Setting build type to '${default_build_type}' as none was specified.")
    set(CMAKE_BUILD_TYPE "${default_build_type}" CACHE
            STRING "Choose the type of build." FORCE)
    # Set the possible values of build type for cmake-gui
    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS
            "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif()

project(LibGeodesic)

set(CMAKE_CXX_STANDARD 17)

find_package(MKL REQUIRED)
find_package(TBB REQUIRED)
set(CGAL_DISABLE_GMP TRUE)

OPTION(DEFINE_MACRO "Option description" ON) # Enabled by default

IF(DEFINE_MACRO)
ADD_DEFINITIONS(-DBUILD_DLL_EXPORTS)
ENDIF(DEFINE_MACRO)

set(BUILD_DLL_EXPORTS)
find_package(CGAL REQUIRED)

include(${CGAL_USE_FILE})

include_directories(${MKL_ROOT}/include ./include ./NurbLibSt ${TBB_ROOT}/include)
link_directories(${MKL_ROOT}/lib/intel64)

find_package(OpenMP)
add_library(LibGeodesic SHARED dllmain.cpp CPointCloud.cpp cmathutilities.cpp main.cpp)
target_include_directories(LibGeodesic PUBLIC TBB::tbb MKL)
target_link_libraries(LibGeodesic PUBLIC mkl_intel_lp64 mkl_tbb_thread mkl_core TBB::tbb)

add_executable(testGeodesic main.cpp)
target_link_libraries(testGeodesic LibGeodesic)

if(OpenMP_CXX_FOUND)
    target_link_libraries(LibGeodesic PUBLIC OpenMP::OpenMP_CXX)
endif()

