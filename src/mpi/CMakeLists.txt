if (OPT_MPI)
    target_sources(lmcore
        PRIVATE
            "${CMAKE_CURRENT_LIST_DIR}/lm_mpi.cpp"
        PUBLIC
            "${CMAKE_CURRENT_LIST_DIR}/lm_mpi.h")
endif()
