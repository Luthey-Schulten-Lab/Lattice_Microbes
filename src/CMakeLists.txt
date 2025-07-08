cmake_minimum_required(VERSION 3.10 FATAL_ERROR)
project(lm VERSION 2.5.0 LANGUAGES C CXX DESCRIPTION "Lattice Microbes")

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if (NOT CMAKE_BUILD_TYPE)
      set(CMAKE_BUILD_TYPE "Release" CACHE STRING
         "Choose the type of build, options are: Debug Release RelWithDebInfo MinSizeRel."
         FORCE)
endif ()

if (UNIX AND NOT APPLE)
    set(LINUX ON)
    set(MACOSX OFF)
elseif (APPLE)
    set(LINUX OFF)
    set(MACOSX ON)
else ()
    message(FATAL_ERROR "Unsupported OS")
endif ()

if (${CMAKE_SYSTEM_PROCESSOR} MATCHES "arm") 
    set(ARM ON)
else ()
    set(ARM OFF)
endif ()

include(GNUInstallDirs)

string(TIMESTAMP COPYRIGHT_DATE "2008-%Y")
string(TIMESTAMP BUILD_TIMESTAMP "%Y-%m-%d %H:%M:%s")
string(TIMESTAMP BUILD_USER "$ENV{USER}")
cmake_host_system_information(RESULT BUILD_HOST QUERY HOSTNAME)

find_package(Git)

if (GIT_FOUND)
    execute_process(
        COMMAND git -P rev-parse --abbrev-ref HEAD
        WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
        OUTPUT_VARIABLE GIT_BRANCH
        RESULT_VARIABLE GIT_NOT_AVAILABLE
        ERROR_QUIET
        OUTPUT_STRIP_TRAILING_WHITESPACE)

    execute_process(
        COMMAND git -P log -1 --format=%h
        WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
        OUTPUT_VARIABLE GIT_COMMIT_HASH
        ERROR_QUIET
        OUTPUT_STRIP_TRAILING_WHITESPACE)
endif ()

if (NOT GIT_NOT_AVAILABLE)
    set(BUILD_ID "${PROJECT_DESCRIPTION} ${PROJECT_VERSION} (git:${GIT_BRANCH}:${GIT_COMMIT_HASH})")
else()
    set(BUILD_ID "${PROJECT_DESCRIPTION} ${PROJECT_VERSION}")
endif()

include(cmake/options.cmake)

include(cmake/python_env.cmake)

find_package(HDF5 REQUIRED COMPONENTS CXX HL)
find_package(Protobuf REQUIRED)
find_package(Threads REQUIRED)

add_library(lmcore "")
set_property(TARGET lmcore PROPERTY POSITION_INDEPENDENT_CODE ON)
target_include_directories(lmcore PUBLIC "${CMAKE_BINARY_DIR}" "${CMAKE_SOURCE_DIR}" "${HDF5_INCLUDE_DIRS}")
target_link_libraries(lmcore PUBLIC Threads::Threads protobuf::libprotobuf ${HDF5_HL_LIBRARIES} ${HDF5_LIBRARIES})
if (CMAKE_BUILD_TYPE MATCHES DEBUG OR CMAKE_BUILD_TYPE MATCHES RelWithDebInfo)
    target_compile_options(lmcore PRIVATE $<$<COMPILE_LANGUAGE:CUDA>:-G -g>)
endif()

if (OPT_PYTHON)
    include(cmake/python.cmake)
    target_include_directories(lmcore PUBLIC ${PYTHON_INCLUDE_DIRS})
    target_link_libraries(lmcore PUBLIC ${PYTHON_LIBRARIES})
endif()


if (OPT_CUDA)
    set(CUDA_ARCHITECTURES "70;75;80;86" CACHE STRING "CUDA architectures to compile for")
    include(cmake/cuda.cmake)
    
    message(STATUS "Final CUDA_ARCHITECTURES: ${CMAKE_CUDA_ARCHITECTURES}")
    
    # include(cmake/cuda.cmake)
    if (PROF_USE_NVTX)
        target_include_directories(lmcore PUBLIC ${NVTX_INCLUDE_DIR})
        target_link_libraries(lmcore PUBLIC ${NVTX_LIBRARY})
    endif ()
    set_property(TARGET lmcore PROPERTY CUDA_SEPARABLE_COMPILATION ON)
    set_property(TARGET lmcore PROPERTY CUDA_STANDARD 11)
endif ()

if (OPT_MPI)
    find_package(MPI REQUIRED)
    target_link_libraries(lmcore PUBLIC MPI::MPI_CXX)
    if (${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION} LESS_EQUAL 3.11)
        message(WARNING "MPI/CUDA builds are currently broken with CMake 3.11 or older, see https://gitlab.kitware.com/cmake/cmake/issues/17929")
    endif()
endif ()


include(builder/CMakeLists.txt)
include(mpi/CMakeLists.txt)
include(builder/CMakeLists.txt)
include(swig/CMakeLists.txt)
include(me/CMakeLists.txt)
include(thread/CMakeLists.txt)
include(core/CMakeLists.txt)
include(reaction/CMakeLists.txt)
include(protobuf/CMakeLists.txt)
include(io/CMakeLists.txt)
include(cme/CMakeLists.txt)
include(rng/CMakeLists.txt)
include(cuda/CMakeLists.txt)
include(rdme/CMakeLists.txt)
include(cmd/CMakeLists.txt)
include(pylm/CMakeLists.txt)
include(jlm/CMakeLists.txt)
include(vmd/CMakeLists.txt)

string(TIMESTAMP BUILD_TIMESTAMP "%Y-%m-%d %H:%M:%S")
string(TIMESTAMP BUILD_USER "$ENV{USER}")
cmake_host_system_information(RESULT BUILD_HOST QUERY HOSTNAME)

#set(BUILD_CONFIG_LIST "${BUILD_ID}")
#list(APPEND BUILD_CONFIG_LIST "Built by ${BUILD_USER}@${BUILD_HOST} at ${BUILD_TIMESTAMP}")
set(BUILD_CONFIG_LIST "C++ compiler: ${CMAKE_CXX_COMPILER_ID} ${CMAKE_CXX_COMPILER_VERSION}")
list(APPEND BUILD_CONFIG_LIST "HDF5 version: ${HDF5_VERSION}")

if (OPT_CUDA)
    list(APPEND BUILD_CONFIG_LIST "CUDA Toolkit Version: ${CMAKE_CUDA_COMPILER_VERSION}")
    list(APPEND BUILD_CONFIG_LIST "CUDA Target Architecture: ${CUDA_COMPILED_ARCH}")
else()
    list(APPEND BUILD_CONFIG_LIST "CUDA support disabled")
endif()

if (OPT_PYTHON)
    list(APPEND BUILD_CONFIG_LIST "Python version: ${PYTHONLIBS_VERSION_STRING}")
    list(APPEND BUILD_CONFIG_LIST "NumPy version: ${_numpy_version}")
    list(APPEND BUILD_CONFIG_LIST "Swig version: ${SWIG_VERSION}")
else()
    list(APPEND BUILD_CONFIG_LIST "Python support disabled")
endif()

if (OPT_MPI)
    list(APPEND BUILD_CONFIG_LIST "MPI version: ${MPI_VERSION}")
else()
    list(APPEND BUILD_CONFIG_LIST "MPI support disabled")
endif()

function(cfg_concat)
    set(BUILD_CONFIGURATION "${BUILD_CONFIGURATION}${ARGN}" PARENT_SCOPE)
endfunction()

cfg_concat("${BUILD_ID}")
cfg_concat("\\nBuilt by ${BUILD_USER}@${BUILD_HOST} at ${BUILD_TIMESTAMP}")

message("============================== Build configuration =============================")
cfg_concat("\\nCompiled features:")
message("Compiled features:")
foreach (var ${BUILD_CONFIG_LIST})
    cfg_concat("\\n * ${var}")
    message(" * ${var}")
endforeach()

cfg_concat("\\nCompiled options:")
foreach (var ${COMPILED_OPTIONS_LIST})
    cfg_concat("\\n * ${var}")
endforeach()

if (COMPILED_OPTIONS_LIST_NONDEFAULT)
    message("Compiled options:")
    foreach (var ${COMPILED_OPTIONS_LIST_NONDEFAULT})
        message(" * ${var}")
    endforeach()
endif()

message("================================================================================")

configure_file(config.h.in config.h)
