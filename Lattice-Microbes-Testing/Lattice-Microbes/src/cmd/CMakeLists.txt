add_library(lmmain STATIC "${CMAKE_CURRENT_LIST_DIR}/common.cpp") 
target_include_directories(lmmain PRIVATE lmcore)
target_link_libraries(lmmain PUBLIC lmcore)

add_executable(lm "${CMAKE_CURRENT_LIST_DIR}/MainSA.cpp") 
target_link_libraries(lm PRIVATE lmmain)

add_executable(lm_setdm "${CMAKE_CURRENT_LIST_DIR}/lm_setdm.cpp")
target_link_libraries(lm_setdm PUBLIC lmcore)

add_executable(lm_setp "${CMAKE_CURRENT_LIST_DIR}/lm_setp.cpp")
target_link_libraries(lm_setp PUBLIC lmcore)

add_executable(lm_setrm "${CMAKE_CURRENT_LIST_DIR}/lm_setrm.cpp")
target_link_libraries(lm_setrm PUBLIC lmcore)

install(TARGETS lm lm_setdm lm_setp lm_setrm RUNTIME DESTINATION bin)

if (OPT_MPI)
     add_executable(mpilm "${CMAKE_CURRENT_LIST_DIR}/MainMPI.cpp")
     target_link_libraries(mpilm PUBLIC lmcore)
     install(TARGETS mpilm RUNTIME DESTINATION bin)
endif ()

if (OPT_PYTHON)
     add_executable(lm_python "${CMAKE_CURRENT_LIST_DIR}/lm_python.cpp" "${CMAKE_BINARY_DIR}/lm_module_pack.h")
     target_include_directories(lm_python PRIVATE ${SWIG_MODULE_pythonlmstatic_REAL_NAME})
     target_link_libraries(lm_python PRIVATE ${SWIG_MODULE_pythonlmstatic_REAL_NAME})
     install(TARGETS lm_python RUNTIME DESTINATION bin)
endif ()
