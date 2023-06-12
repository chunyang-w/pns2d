add_library(CMatrix CMatrix.cpp)
set_target_properties(CMatrix
    PROPERTIES
    CXX_STANDARD 17
    ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib"
)

add_library(helper helper.cpp)
set_target_properties(CMatrix
    PROPERTIES
    CXX_STANDARD 17
    ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib"
)

add_executable(NS_Solver main.cpp)
target_link_libraries(NS_Solver PUBLIC ${MPI_LIBRARIES} CMatrix helper)

set_target_properties(NS_Solver
    PROPERTIES
    CXX_STANDARD 17
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin"
)