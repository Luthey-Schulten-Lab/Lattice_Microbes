if(OPT_PYTHON)
    find_program(CYTHON_EXECUTABLE NAMES cython PATH_SUFFIXES bin DOC "Path to Cython")
    if (CYTHON_EXECUTABLE)
        file(GLOB_RECURSE JLM_FILES LIST_DIRECTORIES FALSE CONFIGURE_DEPENDS "${CMAKE_CURRENT_LIST_DIR}/*")
        add_custom_command(
            OUTPUT "${CMAKE_BINARY_DIR}/jlm.dep"
            DEPENDS "${JLM_FILES}"
            COMMAND "${CMAKE_COMMAND}" -E copy_directory
                    "${CMAKE_SOURCE_DIR}/jlm"
                    "${CMAKE_BINARY_DIR}/jlm"
            COMMAND "${CMAKE_COMMAND}" -E remove "${CMAKE_BINARY_DIR}/jlm/setup.py.in"
            COMMAND "${CMAKE_COMMAND}" -E touch "${CMAKE_BINARY_DIR}/jlm.dep")

        add_custom_target(jlm ALL DEPENDS "${CMAKE_BINARY_DIR}/jlm.dep")

        configure_file("${CMAKE_CURRENT_LIST_DIR}/setup.py.in" "${CMAKE_BINARY_DIR}/jlm/setup.py")

        install(CODE "execute_process(COMMAND \"${PYTHON_EXECUTABLE}\" setup.py install WORKING_DIRECTORY \"${CMAKE_BINARY_DIR}/jlm\")")
    else()
        message(FATAL_ERROR "Cython is required for Python support.")
    endif()
endif()
