if (OPT_VMD)
    find_program(VMD_EXECUTABLE NAMES vmd PATH_SUFFIXES bin DOC "Path to VMD start script")
    if(VMD_EXECUTABLE)
        file(STRINGS "${VMD_EXECUTABLE}" _VMD_ROOT_DIR REGEX "^[^#]*defaultvmddir[ \t]*=")
        string(REGEX REPLACE "(^.*=[^\"]*\"?|\"$)" "" _VMD_ROOT_DIR "${_VMD_ROOT_DIR}")
    else()
        set(_VMD_ROOT_DIR "/usr/local/lib/vmd")
    endif()

    set(VMD_ROOT_DIR "${_VMD_ROOT_DIR}" CACHE PATH "Path to VMD runtime directory") 

    find_path(VMD_PLUGIN_DIR "pdbplugin.so"
                PATHS 
                    "${VMD_ROOT_DIR}/plugins/LINUXAMD64/molfile" 
                    "${VMD_ROOT_DIR}/plugins/LINUX/molfile" 
                    "${VMD_ROOT_DIR}/plugins/MACOSX86/molfile" 
                    "${VMD_ROOT_DIR}/plugins/MACOSX86_64/molfile" 
                DOC
                    "Path to VMD molfile plugin directory"
                NO_DEFAULT_PATH)

    find_path(VMD_PLUGIN_INCLUDE_DIR "vmdplugin.h"
                PATHS 
                    "${VMD_ROOT_DIR}/plugins/include"
                DOC
                    "Path to VMD plugin header directory"
                NO_DEFAULT_PATH)

    mark_as_advanced(VMD_ROOT_DIR VMD_EXECUTABLE)

    if(VMD_PLUGIN_INCLUDE_DIR)
        message(STATUS "Found VMD molfile directory: `${VMD_PLUGIN_DIR}'")
        add_library(lmplugin SHARED "${CMAKE_CURRENT_LIST_DIR}/lmplugin.cpp")
        target_include_directories(lmplugin PUBLIC "${HDF5_INCLUDE_DIRS}" "${VMD_PLUGIN_INCLUDE_DIR}")
        target_link_libraries(lmplugin PUBLIC "${HDF5_HL_LIBRARIES}" "${HDF5_LIBRARIES}")
        set_target_properties(lmplugin
            PROPERTIES 
                PREFIX ""
                INSTALL_RPATH "${CMAKE_PREFIX_PATH}/lib")
        install(TARGETS lmplugin DESTINATION "${VMD_PLUGIN_DIR}")
    else()
        message(WARNING "VMD not found. Please set \"VMD_PLUGIN_INCLUDE_DIR\" and \"VMD_PLUGIN_DIR\".")
    endif()
endif()
