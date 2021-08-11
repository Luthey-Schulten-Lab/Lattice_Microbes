if(OPT_PYTHON)
    file(GLOB_RECURSE PYLM_FILES LIST_DIRECTORIES FALSE CONFIGURE_DEPENDS "${CMAKE_CURRENT_LIST_DIR}/*")
    add_custom_command(
        OUTPUT "${CMAKE_BINARY_DIR}/pylm.dep"
        DEPENDS "${PYLM_FILES}"
        COMMAND "${CMAKE_COMMAND}" -E copy_directory
                "${CMAKE_SOURCE_DIR}/pylm"
                "${CMAKE_BINARY_DIR}/pylm"
        COMMAND "${CMAKE_COMMAND}" -E remove "${CMAKE_BINARY_DIR}/pylm/setup.py.in"
        COMMAND "${CMAKE_COMMAND}" -E touch "${CMAKE_BINARY_DIR}/pylm.dep")

    add_custom_target(pylm ALL DEPENDS "${CMAKE_BINARY_DIR}/pylm.dep")

    configure_file("${CMAKE_CURRENT_LIST_DIR}/setup.py.in" "${CMAKE_BINARY_DIR}/pylm/setup.py")

    install(CODE "execute_process(COMMAND \"${PYTHON_EXECUTABLE}\" setup.py install WORKING_DIRECTORY \"${CMAKE_BINARY_DIR}/pylm\")")
endif()
