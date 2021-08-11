if (OPT_PYTHON)
    set_property(SOURCE "${CMAKE_CURRENT_LIST_DIR}/lm.i"
                 PROPERTY SWIG_MODULE_NAME "lm")

    set_property(SOURCE "${CMAKE_CURRENT_LIST_DIR}/lm.i"
                 PROPERTY CPLUSPLUS ON)

    swig_add_library(pythonlm
         TYPE SHARED
         LANGUAGE python
         SOURCES "${CMAKE_CURRENT_LIST_DIR}/lm.i")

    target_include_directories(${SWIG_MODULE_pythonlm_REAL_NAME}
         PRIVATE
            lmcore
            ${NUMPY_INCLUDE_DIR}
            ${PYTHON_INCLUDE_DIRS})

    target_link_libraries(${SWIG_MODULE_pythonlm_REAL_NAME}
         PRIVATE
            lmcore
            ${PYTHON_LIBRARIES})

    swig_add_library(pythonlmstatic
         TYPE STATIC
         LANGUAGE python
         SOURCES "${CMAKE_CURRENT_LIST_DIR}/lm.i")

    target_include_directories(${SWIG_MODULE_pythonlmstatic_REAL_NAME}
        PUBLIC
            lmcore
            ${NUMPY_INCLUDE_DIR}
            ${PYTHON_INCLUDE_DIRS})

    target_link_libraries(${SWIG_MODULE_pythonlmstatic_REAL_NAME}
        PUBLIC
            lmcore
            ${PYTHON_LIBRARIES})

    add_custom_command(
        OUTPUT "${CMAKE_BINARY_DIR}/lm_module_pack.h"
        DEPENDS ${SWIG_MODULE_pythonlmstatic_REAL_NAME}
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
        COMMAND
            ${PYTHON_EXECUTABLE}
                ${CMAKE_CURRENT_LIST_DIR}/embed_lm_dot_py.py
                ${CMAKE_BINARY_DIR}/lm.py
                ${CMAKE_BINARY_DIR}/lm_module_pack.h)

    add_custom_target(lm_module_pack
        DEPENDS "${CMAKE_BINARY_DIR}/lm_module_pack.h")

    set(SWIG_LM_PYTHON_WRAPPER "${CMAKE_BINARY_DIR}/lm.py")
    set(SWIG_LM_MODULE_OLD_NAME "${CMAKE_BINARY_DIR}/${SWIG_MODULE_pythonlm_REAL_NAME}${CMAKE_SHARED_LIBRARY_SUFFIX}")

    execute_process(
        COMMAND ${PYTHON_EXECUTABLE} "-c"
        "import sysconfig; print(\"_lm\"+sysconfig.get_config_var(\"EXT_SUFFIX\"))"
        OUTPUT_VARIABLE SWIG_LM_MODULE_NAME
        OUTPUT_STRIP_TRAILING_WHITESPACE
    )

    if(APPLE)
        add_custom_target(
            "${SWIG_LM_MODULE_NAME}" ALL
            DEPENDS "${SWIG_MODULE_pythonlm_REAL_NAME}"
            WORKING_DIRECTORY "${CMAKE_BINARY_DIR}"
            COMMAND
                ${CMAKE_COMMAND} -E rename "${SWIG_LM_MODULE_OLD_NAME}" "${SWIG_LM_MODULE_NAME}"
            COMMAND # Not sure if this is necessary for Darwin
                ${CMAKE_INSTALL_NAME_TOOL} -id "@rpath/${SWIG_LM_MODULE_NAME}" "${SWIG_LM_MODULE_NAME}"
            )
    else()
        add_custom_target(
            "${SWIG_LM_MODULE_NAME}" ALL
            DEPENDS "${SWIG_MODULE_pythonlm_REAL_NAME}"
            WORKING_DIRECTORY "${CMAKE_BINARY_DIR}"
            COMMAND
                ${CMAKE_COMMAND} -E rename "${SWIG_LM_MODULE_OLD_NAME}" "${SWIG_LM_MODULE_NAME}"
        )
    endif()

    configure_file(${CMAKE_CURRENT_LIST_DIR}/setup.py.in setup.py)

    install(
        CODE "execute_process(COMMAND ${PYTHON_EXECUTABLE} setup.py install)"
        WORKING_DIRECTORY "${CMAKE_BINARY_DIR}")

endif ()
