include(CheckLanguage)
check_language(CUDA)

if (CMAKE_CUDA_COMPILER)
    set(ENV{CUDAHOSTCXX} "${CMAKE_CXX_COMPILER}")
    enable_language(CUDA)

    if(NOT CUDA_GPU_DETECT_OUTPUT) # taken from FindCUDA/select_compute_arch.cmake
        if(DEFINED ENV{CUDA_ARCHITECTURES})
            message(STATUS "Using CUDA architectures from environment: $ENV{CUDA_ARCHITECTURES}")
            set(DETECTED_CUDA_ARCH $ENV{CUDA_ARCHITECTURES})
        elseif(DEFINED CUDA_ARCHITECTURES)
            message(STATUS "Using manual input CUDA architectures: ${CUDA_ARCHITECTURES}")
            set(DETECTED_CUDA_ARCH ${CUDA_ARCHITECTURES})
        else()
            message(STATUS "Using CUDA architectures from auto-detection")
            set(file ${PROJECT_BINARY_DIR}/detect_cuda_compute_capabilities.cu)

            file(WRITE ${file} "" 
                "#include <cuda_runtime.h>\n"
                "#include <cstdio>\n"
                "int main()\n"
                "{\n"
                "  int count = 0;\n"
                "  if (cudaSuccess != cudaGetDeviceCount(&count)) return -1;\n"
                "  if (count == 0) return -1;\n"
                "  for (int device = 0; device < count; ++device) {\n"
                "    cudaDeviceProp prop;\n"
                "    if (cudaSuccess == cudaGetDeviceProperties(&prop, device)) {\n"
                "      std::printf(\"%d%d\", prop.major, prop.minor);\n"
                "      if (device<count-1) {\n"
                "        std::printf(\";\");\n"
                "      }\n"
                "    }\n"
                "  }\n"
                "  return 0;\n"
                "}\n")

            try_run(run_result compile_result ${PROJECT_BINARY_DIR} ${file}
                    CMAKE_FLAGS "-DINCLUDE_DIRECTORIES=${CUDA_INCLUDE_DIRS}"
                    LINK_LIBRARIES ${CUDA_LIBRARIES}
                    RUN_OUTPUT_VARIABLE compute_capabilities)

            if (run_result EQUAL 0)
                set(CUDA_GPU_DETECT_OUTPUT ${compute_capabilities}
                    CACHE INTERNAL "Returned GPU architectures from detect_gpus tool" FORCE)
                set(DETECTED_CUDA_ARCH ${CUDA_GPU_DETECT_OUTPUT})
            else ()
                message(STATUS "Automatic GPU detection failed. Building for common architectures.")
                set(DETECTED_CUDA_ARCH 70 75 80 86)
                # 90 for H100 Grace Hopper, 86 for A100, 86 for RTX A5000, 75 for RTX 5000
            endif ()
        endif()
    endif()

    set(CUDA_ARCH_LIST "${DETECTED_CUDA_ARCH}" CACHE STRING "CUDA arch to compile code for. Specify a semicolon delimited list, e.g. 35;52")
    message(STATUS "CUDA_ARCH_LIST: ${CUDA_ARCH_LIST}")
    list(SORT CUDA_ARCH_LIST)
    message(STATUS "CUDA_ARCH_LIST after sorting: ${CUDA_ARCH_LIST}")
    set(CUDA_ARCH_FLAGS "" CACHE STRING "CUDA arch command line options")
    mark_as_advanced(CUDA_ARCH_FLAGS)

    foreach(x ${CUDA_ARCH_LIST})
        message(STATUS "Adding CUDA architecture: compute_${x}, sm_${x}")
        set(CUDA_ARCH_FLAGS "${CUDA_ARCH_FLAGS} -gencode arch=compute_${x},code=sm_${x}")
        if(CUDA_COMPILED_ARCH)
            set(CUDA_COMPILED_ARCH "${CUDA_COMPILED_ARCH}, sm_${x}")
        else()
            set(CUDA_COMPILED_ARCH "sm_${x}")
        endif()
    endforeach()


    if(CUDA_COMPILED_ARCH)
        message(STATUS "Building for CUDA architectures: ${CUDA_COMPILED_ARCH}")
    endif()
    set(CMAKE_CUDA_FLAGS "${CUDA_ARCH_FLAGS}")

    if (PROF_USE_NVTX)
        find_library(NVTX_LIBRARY nvToolsExt HINT ${CMAKE_CUDA_IMPLICIT_LINK_DIRECTORIES})
        set(NVTX_INCLUDE_DIR "${CMAKE_CUDA_TOOLKIT_INCLUDE_DIRECTORIES}")
    endif()
else()
    message(WARNING "NVCC not found, disabling CUDA support")
    set(OPT_CUDA OFF)
endif()
