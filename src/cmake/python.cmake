find_package(PythonInterp 3.5)

function(verify_python_module name module)
    if (OPT_PYTHON)
        execute_process(
            COMMAND ${PYTHON_EXECUTABLE} "-c" "import ${module}"
            RESULT_VARIABLE _${module}_status
            ERROR_QUIET
            OUTPUT_STRIP_TRAILING_WHITESPACE)

        if(_${module}_status)
            message(WARNING "${name} is required to build with Python support. Python disabled.")
            set(OPT_PYTHON OFF PARENT_SCOPE)
        else()
            message(STATUS "Found ${name}")
        endif()
    endif()
endfunction()

if (PythonInterp_FOUND)
    find_package(PythonLibs ${PYTHON_VERSION_MAJOR}.${PYTHON_VERSION_MINOR} REQUIRED)
    find_package(SWIG REQUIRED)
    include(UseSWIG)

    execute_process(
        COMMAND ${PYTHON_EXECUTABLE} "-c" "import numpy; print(numpy.get_include())"
        RESULT_VARIABLE _numpy_status
        OUTPUT_VARIABLE _numpy_location
        ERROR_QUIET
        OUTPUT_STRIP_TRAILING_WHITESPACE)

    execute_process(
        COMMAND ${PYTHON_EXECUTABLE} "-c" "import numpy; print(numpy.__version__)"
        RESULT_VARIABLE _numpy_status
        OUTPUT_VARIABLE _numpy_version
        ERROR_QUIET
        OUTPUT_STRIP_TRAILING_WHITESPACE)

    if(NOT _numpy_status)
        set(NUMPY_INCLUDE_DIR ${_numpy_location} CACHE STRING "Location of NumPy")
    else()
        message(WARNING "NumPy is required to build with Python support. Python disabled.")
        set(OPT_PYTHON OFF)
    endif()

    include(FindPackageHandleStandardArgs)
    find_package_handle_standard_args(NumPy
        FOUND_VAR NumPy_FOUND
        REQUIRED_VARS NUMPY_INCLUDE_DIR
        VERSION_VAR _numpy_version)

    verify_python_module("Cython" cython)
    verify_python_module("Setuptools" setuptools)
else()
    message(WARNING "Python 3.5 or newer not found. Python support disabled")
    set(OPT_PYTHON OFF)
endif()
