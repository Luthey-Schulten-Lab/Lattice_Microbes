target_sources(lmcore
    PRIVATE
        "${CMAKE_CURRENT_LIST_DIR}/ByteLattice.cpp"
        "${CMAKE_CURRENT_LIST_DIR}/IntLattice.cpp"
        "${CMAKE_CURRENT_LIST_DIR}/Lattice.cpp"
        "${CMAKE_CURRENT_LIST_DIR}/NextSubvolumeSolver.cpp"
        "${CMAKE_CURRENT_LIST_DIR}/RDMESolver.cpp"
    PUBLIC
        "${CMAKE_CURRENT_LIST_DIR}/NextSubvolumeSolver.h"
        "${CMAKE_CURRENT_LIST_DIR}/RDMESolver.h"
        "${CMAKE_CURRENT_LIST_DIR}/ByteLattice.h"
        "${CMAKE_CURRENT_LIST_DIR}/Lattice.h"
        "${CMAKE_CURRENT_LIST_DIR}/IntLattice.h")

if (OPT_CUDA)
    target_sources(lmcore
        PRIVATE
            "${CMAKE_CURRENT_LIST_DIR}/CudaByteLattice.cu"
            "${CMAKE_CURRENT_LIST_DIR}/CudaIntLattice.cu"
            "${CMAKE_CURRENT_LIST_DIR}/IntMpdRdmeSolver.cu"
            "${CMAKE_CURRENT_LIST_DIR}/MGPUMpdRdmeSolver.cu"
            "${CMAKE_CURRENT_LIST_DIR}/MpdRdmeSolver.cu"
            "${CMAKE_CURRENT_LIST_DIR}/GPUMapper/MultiGPUMapper.cu"
            "${CMAKE_CURRENT_LIST_DIR}/GPUMapper/ZDivMultiGPUMapper.cu"
            "${CMAKE_CURRENT_LIST_DIR}/GPUMapper/osx_barrier.cu"
        PUBLIC
            "${CMAKE_CURRENT_LIST_DIR}/CudaByteLattice.h"
            "${CMAKE_CURRENT_LIST_DIR}/CudaIntLattice.h"
            "${CMAKE_CURRENT_LIST_DIR}/IntMpdRdmeSolver.h"
            "${CMAKE_CURRENT_LIST_DIR}/MGPUMpdRdmeSolver.h"
            "${CMAKE_CURRENT_LIST_DIR}/MpdRdmeSolver.h")

    if (OPT_MPI)
        target_sources(lmcore
            PRIVATE
                "${CMAKE_CURRENT_LIST_DIR}/MPIMpdRdmeSolver.cu"
            PUBLIC
                "${CMAKE_CURRENT_LIST_DIR}/MPIMpdRdmeSolver.h")
    endif ()
endif()
