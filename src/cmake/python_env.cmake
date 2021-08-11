find_program(PYTHON_EXECUTABLE python PATHS ENV PATH NO_DEFAULT_PATH) # force CMake to use the venv binary if it exists
if(OPT_PYTHON)
    execute_process(
        COMMAND ${PYTHON_EXECUTABLE} "-c" "import sys; print(sys.exec_prefix)"
        RESULT_VARIABLE _python_prefix_status
        OUTPUT_VARIABLE _python_prefix
        ERROR_QUIET
        OUTPUT_STRIP_TRAILING_WHITESPACE)
    set(DETECTED_PREFIX "${_python_prefix}"
        CACHE INTERNAL "Prefix inferred from python binary" FORCE)
    if(NOT _python_prefix_status)
        set(CMAKE_INSTALL_PREFIX "${DETECTED_PREFIX}" CACHE PATH "Install prefix" FORCE)
        set(CMAKE_PREFIX_PATH "${DETECTED_PREFIX}" CACHE PATH "Virtual environment libraries" FORCE)
        message(STATUS "Inferring install prefix from '${DETECTED_PREFIX}'")
    endif()

    # If we're compiling using anaconda gcc, if we don't set the sysroot g++ will try to link the wrong librt
    if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
        execute_process(
            COMMAND ${CMAKE_C_COMPILER} "-print-sysroot"
            RESULT_VARIABLE _sysroot_status
            OUTPUT_VARIABLE _sysroot_prefix
            ERROR_QUIET
            OUTPUT_STRIP_TRAILING_WHITESPACE)
        if(NOT _sysroot_status AND NOT "${_sysroot_prefix}" STREQUAL "")
            set(DETECTED_SYSROOT "${_sysroot_prefix}"
                CACHE INTERNAL "Inferred GCC sysroot" FORCE)
            set(CMAKE_SYSROOT "${DETECTED_SYSROOT}" CACHE PATH "GCC Sysroot" FORCE)
            message(STATUS "Cross-compiling: GCC sysroot is '${DETECTED_SYSROOT}'")
        endif()
    endif()
endif()
