# Gurobi ships headers and libraries without a CMake package configuration.
set(GUROBI_HOME "$ENV{GUROBI_HOME}" CACHE PATH "Gurobi platform directory (containing include and lib)")
find_path(GUROBI_INCLUDE_DIR gurobi_c++.h HINTS "${GUROBI_HOME}" PATH_SUFFIXES include)
if(GUROBI_INCLUDE_DIR)
    file(STRINGS "${GUROBI_INCLUDE_DIR}/gurobi_c.h" _gurobi_version
        REGEX "^#define GRB_VERSION_(MAJOR|MINOR) +[0-9]+")
    string(REGEX MATCH "GRB_VERSION_MAJOR +([0-9]+)" _match "${_gurobi_version}")
    set(_gurobi_major "${CMAKE_MATCH_1}")
    string(REGEX MATCH "GRB_VERSION_MINOR +([0-9]+)" _match "${_gurobi_version}")
    set(_gurobi_suffix "${_gurobi_major}${CMAKE_MATCH_1}")
endif()
find_library(GUROBI_LIBRARY NAMES gurobi${_gurobi_suffix}
    HINTS "${GUROBI_HOME}" PATH_SUFFIXES lib)
if(WIN32)
    if(NOT MSVC)
        message(FATAL_ERROR "The installed Windows Gurobi C++ libraries require MSVC. Use build.ps1 or a Visual Studio CMake generator.")
    endif()
    if(NOT CMAKE_SIZEOF_VOID_P EQUAL 8)
        message(FATAL_ERROR "Gurobi win64 requires an x64 build (use -A x64).")
    endif()
    find_library(GUROBI_CXX_RELEASE NAMES gurobi_c++md2017 HINTS "${GUROBI_HOME}" PATH_SUFFIXES lib)
    find_library(GUROBI_CXX_DEBUG NAMES gurobi_c++mdd2017 HINTS "${GUROBI_HOME}" PATH_SUFFIXES lib)
    find_file(GUROBI_DLL NAMES gurobi${_gurobi_suffix}.dll HINTS "${GUROBI_HOME}" PATH_SUFFIXES bin)
    set(_gurobi_required GUROBI_CXX_RELEASE GUROBI_CXX_DEBUG GUROBI_DLL)
else()
    find_library(GUROBI_CXX_RELEASE NAMES gurobi_c++ HINTS "${GUROBI_HOME}" PATH_SUFFIXES lib)
    set(GUROBI_CXX_DEBUG "${GUROBI_CXX_RELEASE}")
    set(_gurobi_required GUROBI_CXX_RELEASE)
endif()
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(GUROBI REQUIRED_VARS
    GUROBI_INCLUDE_DIR GUROBI_LIBRARY ${_gurobi_required})
if(GUROBI_FOUND AND NOT TARGET GUROBI::CXX)
    add_library(GUROBI::CXX UNKNOWN IMPORTED)
    set_target_properties(GUROBI::CXX PROPERTIES
        IMPORTED_LOCATION "${GUROBI_CXX_RELEASE}"
        IMPORTED_LOCATION_DEBUG "${GUROBI_CXX_DEBUG}"
        INTERFACE_INCLUDE_DIRECTORIES "${GUROBI_INCLUDE_DIR}"
        INTERFACE_LINK_LIBRARIES "${GUROBI_LIBRARY}")
endif()
