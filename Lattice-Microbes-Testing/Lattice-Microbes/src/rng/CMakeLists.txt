target_sources(lmcore
    PRIVATE
        "${CMAKE_CURRENT_LIST_DIR}/RandomGenerator.cpp"
        "${CMAKE_CURRENT_LIST_DIR}/XORShift.cpp"
    PUBLIC
        "${CMAKE_CURRENT_LIST_DIR}/XORShift.h"
        "${CMAKE_CURRENT_LIST_DIR}/RandomGenerator.h")

if (OPT_CUDA)
    target_sources(lmcore
        PRIVATE
            "${CMAKE_CURRENT_LIST_DIR}/XORWow.cu"
        PUBLIC
            "${CMAKE_CURRENT_LIST_DIR}/XORWow.h")
endif ()
