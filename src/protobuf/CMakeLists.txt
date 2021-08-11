protobuf_generate_cpp(PROTO_SRCS PROTO_HDRS
    "${CMAKE_CURRENT_LIST_DIR}/DiffusionModel.proto"
    "${CMAKE_CURRENT_LIST_DIR}/FirstPassageTimes.proto"
    "${CMAKE_CURRENT_LIST_DIR}/Lattice.proto"
    "${CMAKE_CURRENT_LIST_DIR}/ParameterValues.proto"
    "${CMAKE_CURRENT_LIST_DIR}/ReactionModel.proto"
    "${CMAKE_CURRENT_LIST_DIR}/SimulationParameters.proto"
    "${CMAKE_CURRENT_LIST_DIR}/SpatialModel.proto"
    "${CMAKE_CURRENT_LIST_DIR}/SpeciesCounts.proto"
)

target_sources(lmcore
    PRIVATE
        ${PROTO_SRCS}
    PUBLIC
        ${PROTO_HDRS}
)

