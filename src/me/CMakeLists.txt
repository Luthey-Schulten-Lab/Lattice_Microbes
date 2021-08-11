target_sources(lmcore
    PRIVATE
        "${CMAKE_CURRENT_LIST_DIR}/MESolver.cpp"
        "${CMAKE_CURRENT_LIST_DIR}/MESolverFactory.cpp"
    PUBLIC
        "${CMAKE_CURRENT_LIST_DIR}/MESolverFactory.h"
        "${CMAKE_CURRENT_LIST_DIR}/MESolver.h"
)
