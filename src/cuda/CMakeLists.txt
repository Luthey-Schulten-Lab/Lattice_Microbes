if (OPT_CUDA)
    target_sources(lmcore
        PRIVATE
            "${CMAKE_CURRENT_LIST_DIR}/lm_cuda.cu"
            "${CMAKE_CURRENT_LIST_DIR}/constant.cu"
        PUBLIC
            "${CMAKE_CURRENT_LIST_DIR}/lm_cuda.h"
            "${CMAKE_CURRENT_LIST_DIR}/ldg.h"
            "${CMAKE_CURRENT_LIST_DIR}/constant.cuh")
endif ()
